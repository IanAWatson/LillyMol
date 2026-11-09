"""Run via LillyMol run_python.sh after staging the rebuilt nanobind module."""
import base64
import copy
import importlib.util
import io
import json
from pathlib import Path
import re
import sys
import tempfile
import unittest
from unittest import mock

BIN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BIN))
import activity_report
from dicer_fragment_activity_html import main as html_main
from lillymol import MolFromSmiles


class FragmentHighlightingTest(unittest.TestCase):
    def setUp(self):
        activity_report._fragment_query.cache_clear()

    def test_attachment_sites_match_and_json_is_unchanged(self):
        examples = [dict(smiles='CC(O)N', identifier='parent', activity=3)]
        before = copy.deepcopy(examples)
        labelled = activity_report.labelled_fragment_examples('C[1CH2]O', examples)
        parent = MolFromSmiles(labelled[0]['labelled_smiles'])
        self.assertEqual(parent.number_isotopic_atoms(), 3)
        for atom in range(parent.natoms()):
            self.assertEqual(parent.isotope(atom), 0 if parent.atomic_number(atom) == 7 else 1)
        self.assertEqual(labelled[0]['nhits'], 1)
        self.assertEqual(examples, before)

    def test_union_of_all_embeddings(self):
        examples = [dict(smiles='CC(O)N.CC(O)N', identifier='two', activity=3)]
        labelled = activity_report.labelled_fragment_examples('C[1CH2]O', examples)
        self.assertEqual(MolFromSmiles(labelled[0]['labelled_smiles']).number_isotopic_atoms(), 6)
        self.assertEqual(labelled[0]['nhits'], 2)

    def test_existing_labels_are_preserved(self):
        labelled = dict(smiles='CC(O)N', labelled_smiles='C[1CH](O)N', identifier='marked', activity=1)
        isotopic = dict(smiles='[13CH3]C(O)N', identifier='isotope', activity=2)
        result = activity_report.labelled_fragment_examples('C[1CH2]O', [labelled, isotopic])
        self.assertIs(result[0], labelled)
        self.assertIs(result[1], isotopic)

    def test_no_match_does_not_fall_back_to_unrestricted_search(self):
        examples = [dict(smiles='CCC(O)N', identifier='wrong attachment', activity=1)]
        result = activity_report.labelled_fragment_examples('C[1CH2]O', examples)
        self.assertNotIn('labelled_smiles', result[0])
        self.assertEqual(result[0]['highlight_note'], 'Fragment match unavailable')

    def test_unlabelled_fragment_and_existing_example_skip_bindings(self):
        examples = [dict(smiles='CCO', identifier='plain', activity=1)]
        with mock.patch.object(activity_report, '_fragment_matching_modules', side_effect=AssertionError):
            self.assertIs(activity_report.labelled_fragment_examples('CCO', examples), examples)
            examples[0]['labelled_smiles'] = '[1CH3]CO'
            self.assertIs(activity_report.labelled_fragment_examples('[1CH3]CO', examples), examples)

    def test_query_is_built_once_per_fragment(self):
        examples = [dict(smiles='CC(O)N', identifier='parent', activity=1)]
        activity_report.labelled_fragment_examples('C[1CH2]O', examples)
        activity_report.labelled_fragment_examples('C[1CH2]O', examples)
        self.assertEqual(activity_report._fragment_query.cache_info().misses, 1)
        self.assertEqual(activity_report._fragment_query.cache_info().hits, 1)

    def test_missing_bindings_keeps_examples(self):
        examples = [dict(smiles='CC(O)N', identifier='parent', activity=1)]
        with mock.patch.object(activity_report, '_fragment_matching_modules', return_value=None):
            self.assertIs(activity_report.labelled_fragment_examples('C[1CH2]O', examples), examples)

    def test_old_bindings_report_feature_unavailable(self):
        query = mock.Mock()
        query.build_from_molecule.side_effect = TypeError('old keyword API')
        with mock.patch.object(activity_report, '_fragment_matching_modules',
                               return_value=(MolFromSmiles, lambda: query, object)), \
             mock.patch.object(activity_report, '_warn_attachment_binding_unavailable') as warn:
            self.assertIsNone(activity_report._fragment_query('C[1CH2]O'))
            self.assertIsNone(activity_report._fragment_query('C[1CH2]O'))
            warn.assert_called_once()

    @unittest.skipUnless(importlib.util.find_spec('rdkit'), 'requires RDKit')
    def test_html_command_contains_highlighted_parent(self):
        row = dict(name='fragment', delta_median=2, delta_mean=2,
                   matched=dict(summary=dict(n=1, median=3, mean=3), histogram=[0, 1]),
                   unmatched=dict(summary=dict(n=1, median=1, mean=1), histogram=[1, 0]))
        data = dict(metadata=dict(feature_kind='fragment', queries=1, activity='activity',
                                  molecules_read_with_activity=2),
                    queries=[row])
        data['global'] = dict(summary=dict(median=2, mean=2), histogram=[1, 1])
        row['source'] = 'C[1CH2]O'
        row['examples'] = [dict(smiles='CC(O)N', identifier='parent', activity=3)]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'activity.json'
            path.write_text(json.dumps(data))
            original = path.read_bytes()
            output = io.StringIO()
            with mock.patch('sys.stdout', output):
                self.assertEqual(html_main([str(path), '-o', '-']), 0)
            document = output.getvalue()
            encoded = re.search(r'data:image/svg\+xml;base64,([^\"]+)', document).group(1)
            svg = base64.b64decode(encoded).decode()
            self.assertIn('<ellipse', svg)
            self.assertIn('#33CC66', svg)
            self.assertIn('hits 1', document)
            self.assertEqual(path.read_bytes(), original)


if __name__ == '__main__':
    unittest.main()
