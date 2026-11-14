"""Run with PYTHON=/path/to/matching/python $LILLYMOL_HOME/run_python.sh this_file."""
from contextlib import redirect_stderr
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

BIN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BIN))
import substructure_exploration as exploration
from lillymol import MolFromSmiles, Standardise


class RepresentativeExamplesTest(unittest.TestCase):
    def record(self, identifier, activity, smiles):
        return exploration.MoleculeRecord(smiles, identifier, activity, MolFromSmiles(smiles))

    def test_midpoint_is_activity_range_not_median(self):
        records = [self.record(str(activity), activity, 'C') for activity in [1, 2, 6, 0, 10]]
        chosen = exploration.representative_records(records, 3)
        self.assertEqual([r.activity for r in chosen], [0, 6, 10])

    def test_ties_prefer_smallest_molecule(self):
        records = [self.record('low-large', 0, 'CCC'), self.record('low-small', 0, 'C'),
                   self.record('high-large', 10, 'CCC'), self.record('high-small', 10, 'C'),
                   self.record('middle-large', 4, 'CCC'), self.record('middle-small', 6, 'CC')]
        chosen = exploration.representative_records(records, 3)
        self.assertEqual([r.identifier for r in chosen], ['low-small', 'middle-small', 'high-small'])
        self.assertEqual([r.identifier for r in exploration.representative_records(records[::-1], 3)],
                         [r.identifier for r in chosen])

    def test_constant_activity_and_small_sets(self):
        records = [self.record('large', 5, 'CCC'), self.record('small', 5, 'C'),
                   self.record('medium', 5, 'CC'), self.record('largest', 5, 'CCCC')]
        self.assertEqual([r.identifier for r in exploration.representative_records(records, 3)],
                         ['small', 'medium', 'large'])
        self.assertEqual(exploration.representative_records(records, 0), [])
        self.assertEqual(exploration.representative_records([], 3), [])
        self.assertEqual(len(exploration.representative_records(records[:2], 3)), 2)

    def test_configurable_example_limit(self):
        records = [self.record(str(i), i, 'C') for i in range(5)]
        for limit, expected in [(1, [0]), (2, [0, 4]), (3, [0, 2, 4]),
                                (5, [0, 1, 2, 3, 4]), (10, [0, 1, 2, 3, 4])]:
            with self.subTest(limit=limit):
                self.assertEqual([r.activity for r in exploration.representative_records(records, limit)],
                                 expected)

    def test_analysis_selects_then_labels_examples(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            (work / 'input.smi').write_text(
                'NCCC middle-large\nNCCCC low-large\nNCCC high-large\n'
                'NC middle-small\nN low-small\nN high-small\n')
            (work / 'activity.txt').write_text(
                'ID activity\nmiddle-large 5\nlow-large 0\nhigh-large 10\n'
                'middle-small 5\nlow-small 0\nhigh-small 10\n')
            (work / 'queries.textproto').write_text('query { name: "nitrogen" smarts: "N" }\n')
            output = work / 'result.json'
            self.assertEqual(exploration.main([
                'analyse', '--config', str(work / 'queries.textproto'),
                '-A', str(work / 'activity.txt'), '-o', str(output), str(work / 'input.smi')]), 0)
            row = json.loads(output.read_text())['queries'][0]
            self.assertEqual(row['matched']['summary']['n'], 6)
            self.assertEqual([e['identifier'] for e in row['examples']],
                             ['low-small', 'middle-small', 'high-small'])
            for example in row['examples']:
                self.assertEqual(example['nhits'], 1)
                self.assertEqual(MolFromSmiles(example['labelled_smiles']).number_isotopic_atoms(), 1)


class CatalogueTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.queries = exploration.build_queries(exploration.read_config(
            str(BIN.parent / 'data/substructure_exploration/queries')))
        cls.by_name = {q.name: q.query for q in cls.queries}
        cls.standardiser = Standardise()
        cls.standardiser.activate_all()

    def count(self, name, smiles):
        mol = MolFromSmiles(smiles)
        self.assertIsNotNone(mol, smiles)
        self.standardiser.process(mol)
        matches = self.by_name[name].substructure_search_matches(mol)
        return len(matches) if matches is not None else 0

    def test_catalogue_examples(self):
        examples = json.loads((BIN / 'test/substructure_exploration/examples.json').read_text())
        self.assertEqual([q.name for q in self.queries], [e['name'] for e in examples])
        for example in examples:
            with self.subTest(feature=example['name']):
                self.assertGreater(self.count(example['name'], example['positive']), 0)
                self.assertEqual(self.count(example['name'], example['negative']), 0)

    def test_fluorine_features(self):
        for smiles, cf3, other in [
            ('FC(F)(F)c1ccc(F)cc1', 1, 0), ('CC(F)F', 0, 2),
            ('FC(F)(F)F', 4, 4), ('CC(F)(F)F', 1, 0),
            ('Fc1ccccc1', 0, 0),
        ]:
            with self.subTest(smiles=smiles):
                self.assertEqual(self.count('cf3', smiles), cf3)
                self.assertEqual(self.count('fluorine_not_cf3', smiles), other)
        self.assertEqual(self.count('aryl_fluoride', 'Fc1ccccc1'), 1)
        self.assertEqual(self.count('aryl_fluoride', 'FC(F)(F)c1ccc(F)cc1'), 1)
        for halogen in ('Cl', 'Br', 'I'):
            self.assertEqual(self.count('heavy_halogen', 'CC' + halogen), 1)
        self.assertEqual(self.count('heavy_halogen', 'CCF'), 0)

    def test_amine_exclusions(self):
        for smiles in ('CC(=O)N', 'CS(=O)(=O)N', 'NC(=O)N', 'NC(=N)N', 'CNO', 'CNN'):
            with self.subTest(smiles=smiles):
                self.assertEqual(self.count('primary_amine', smiles), 0)

    def test_amine_environments(self):
        cases = {
            'primary_amine': ('CCN', [
                'NC(=S)C', 'NC#N', 'NP(=O)(O)O', 'NC(=N)C', 'NN', 'NO', 'NSC']),
            'secondary_amine': ('CCNCC', [
                'CNC(=O)C', 'CNC#N', 'CNNC', 'CNOC',
                'CNP(=O)(O)O', 'CNSC', 'CNC(=N)C']),
            'tertiary_amine': ('CN(C)C', [
                'CN(C)C(=O)C', 'CN(C)C#N', 'CN(C)NC', 'CN(C)OC',
                'CN(C)P(=O)(O)O', 'CN(C)SC', 'CN(C)C(=N)C']),
            'cyclic_amine': ('N1CCCCC1', [
                'O=C1NCCC1', 'N1(C#N)CCCC1', 'ON1CCCC1', 'NN1CCCC1',
                'O=P(O)(O)N1CCCC1', 'CSN1CCCC1', 'N=C1NCCC1']),
        }
        for name, (positive, negatives) in cases.items():
            self.assertEqual(self.count(name, positive), 1)
            for smiles in negatives:
                with self.subTest(feature=name, smiles=smiles):
                    self.assertEqual(self.count(name, smiles), 0)
        # Exclusions apply at each nitrogen, not to the entire molecule.
        self.assertEqual(self.count('primary_amine', 'NCCNC(=O)C'), 1)
        self.assertEqual(self.count('secondary_amine', 'NCCNC(=O)C'), 0)
        # A carbonyl two bonds away does not disqualify an amine.
        self.assertEqual(self.count('primary_amine', 'NCC(=O)C'), 1)

    def test_methoxy(self):
        for smiles in ('COc1ccccc1', 'COCC'):
            with self.subTest(smiles=smiles):
                self.assertEqual(self.count('methoxy', smiles), 1)
        for smiles in ('CC(=O)OC', 'CCOc1ccccc1', 'CO', 'COC=C'):
            with self.subTest(smiles=smiles):
                self.assertEqual(self.count('methoxy', smiles), 0)

    def test_restricted_amide_and_ether_definitions(self):
        for smiles, amide in [('CC(=O)NC', 1), ('NC(=O)N', 0),
                              ('COC(=O)N', 0), ('NC=O', 0)]:
            with self.subTest(smiles=smiles):
                self.assertEqual(self.count('amide', smiles), amide)
        for smiles, acyclic, cyclic in [('COC', 1, 0), ('COc1ccccc1', 0, 0),
                                        ('O1CCCC1', 0, 1), ('O=C1OCCC1', 0, 0)]:
            with self.subTest(smiles=smiles):
                self.assertEqual(self.count('ether', smiles), acyclic)
                self.assertEqual(self.count('ether_cyclic', smiles), cyclic)

    def test_standardisation(self):
        self.assertEqual(self.count('acid', 'CC(=O)[O-]'), 1)
        self.assertEqual(self.count('primary_amine', 'CC[NH3+]'), 1)
        self.assertEqual(self.count('quaternary_ammonium', 'C[N+](C)(C)C'), 1)

    def test_minimum_matched_molecules(self):
        with tempfile.TemporaryDirectory(prefix='exploration support ') as tmp:
            work = Path(tmp)
            (work / 'input.smi').write_text('CC first\nCN second\n')
            (work / 'activity.txt').write_text('ID activity\nfirst 1\nsecond 3\n')
            (work / 'queries.textproto').write_text(
                'query { name: "carbon" smarts: "C" }\n'
                'query { name: "nitrogen" smarts: "N" }\n'
                'query { name: "fluorine" smarts: "F" }\n')
            env = dict(os.environ, PYTHON=sys.executable)
            command = [str(BIN / 'substructure_exploration.sh'), 'analyse',
                       '--config', 'queries.textproto', '-A', 'activity.txt',
                       '--max_examples', '0', '-o', 'result.json', '--html', 'report.html', 'input.smi']
            all_names = {'carbon', 'nitrogen', 'fluorine'}
            for threshold, extra, expected in [
                (1, [], {'carbon', 'nitrogen'}),
                (0, ['--min_matches', '0'], all_names),
                (2, ['--min_matches', '2'], {'carbon'}),
                (3, ['--min_matches', '3'], set()),
            ]:
                with self.subTest(threshold=threshold):
                    subprocess.run(command + extra, cwd=work, env=env,
                                   check=True, capture_output=True, text=True)
                    data = json.loads((work / 'result.json').read_text())
                    self.assertEqual({q['name'] for q in data['queries']}, expected)
                    self.assertEqual(data['metadata']['queries'], 3)
                    self.assertEqual(data['metadata']['queries_reported'], len(expected))
                    self.assertEqual(data['metadata']['min_matches'], threshold)
                    self.assertEqual(data['global']['summary']['n'], 2)
                    report = (work / 'report.html').read_text()
                    self.assertIn(f'{len(expected)} of 3 queries', report)
                    for name in all_names:
                        self.assertEqual(f'<strong>{name}</strong>' in report, name in expected)
                    for query in data['queries']:
                        self.assertEqual(query['examples'], [])
                        self.assertEqual(query['matched']['summary']['n'],
                                         {'carbon': 2, 'nitrogen': 1, 'fluorine': 0}[query['name']])
            result = subprocess.run(command + ['--min_matches', '-1'], cwd=work, env=env,
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('non-negative integer', result.stderr)

    def test_wrapper_and_reports(self):
        with tempfile.TemporaryDirectory(prefix='exploration test ') as tmp:
            work = Path(tmp)
            (work / 'input.smi').write_text('CC(=O)[O-] acid\nCCN amine\n')
            (work / 'activity.txt').write_text('ID activity\nacid 1\namine 3\n')
            env = dict(os.environ, PYTHON=sys.executable)
            command = [str(BIN / 'substructure_exploration.sh'), 'analyse',
                       '-A', 'activity.txt', '-o', 'result.json', '--html', 'report.html', 'input.smi']
            subprocess.run(command, cwd=work, env=env, check=True, capture_output=True, text=True)
            data = json.loads((work / 'result.json').read_text())
            self.assertTrue(data['metadata']['standardize'])
            self.assertEqual(data['metadata']['queries'], len(self.queries))
            acid = next(q for q in data['queries'] if q['name'] == 'acid')
            self.assertEqual(acid['matched']['summary']['n'], 1)
            report_html = (work / 'report.html').read_text()
            self.assertIn('<html', report_html)
            if importlib.util.find_spec('rdkit'):
                self.assertIn('<img class="structure"', report_html)
                self.assertIn('src="data:image/svg+xml;base64,', report_html)
                self.assertIn('<summary>SMILES</summary>', report_html)
            subprocess.run(command + ['--no-standardize', '--min_matches', '0'], cwd=work, env=env,
                           check=True, capture_output=True, text=True)
            data = json.loads((work / 'result.json').read_text())
            acid = next(q for q in data['queries'] if q['name'] == 'acid')
            self.assertFalse(data['metadata']['standardize'])
            self.assertEqual(acid['matched']['summary']['n'], 0)
            subprocess.run([str(BIN / 'substructure_exploration.sh'), 'html',
                            'result.json', '-o', 'regenerated.html'], cwd=work,
                           env=env, check=True, capture_output=True, text=True)
            self.assertIn('<html', (work / 'regenerated.html').read_text())
            ranked = [str(BIN / 'substructure_exploration.sh'), 'html', 'result.json',
                      '-o', 'ranked.html', '--sort', 'abs_median', '--top_n', '1',
                      '--min_matches', '1']
            subprocess.run(ranked, cwd=work, env=env, check=True, capture_output=True, text=True)
            self.assertEqual((work / 'ranked.html').read_text().count('<tr data-search='), 1)
            subprocess.run(ranked + ['--min_unmatched', '2'], cwd=work, env=env,
                           check=True, capture_output=True, text=True)
            self.assertEqual((work / 'ranked.html').read_text().count('<tr data-search='), 0)



class QuerySelectionTest(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory(prefix='query selections ')
        self.addCleanup(self._tmp.cleanup)
        self._directory = Path(self._tmp.name)

    def write_query(self, filename, name, smarts):
        path = self._directory / filename
        path.write_text(f'name: "{name}"\nquery {{ smarts: "{smarts}" }}\n')
        return path

    def load(self, path):
        return exploration.build_queries(exploration.read_config(str(path)))

    def assert_bad(self, path, message):
        errors = io.StringIO()
        with redirect_stderr(errors), self.assertRaises(SystemExit) as raised:
            self.load(path)
        self.assertEqual(raised.exception.code, 1)
        self.assertIn(message, errors.getvalue())

    def test_relative_paths_spaces_comments_and_order(self):
        nitrogen = self.write_query('query nitrogen.textproto', 'nitrogen', 'N')
        carbon = self.write_query('carbon.textproto', 'carbon', 'C')
        selection = self._directory / 'queries'
        selection.write_text('# A selected subset\n\nquery nitrogen.textproto\ncarbon.textproto\n')
        queries = self.load(selection)
        self.assertEqual([q.name for q in queries], ['nitrogen', 'carbon'])
        self.assertEqual([q.source for q in queries], [str(nitrogen), str(carbon)])
        molecule = MolFromSmiles('CN')
        self.assertEqual([len(q.query.substructure_search_matches(molecule)) for q in queries], [1, 1])

    def test_environment_variable_paths(self):
        nitrogen = self.write_query('n.textproto', 'nitrogen', 'N')
        carbon = self.write_query('c.textproto', 'carbon', 'C')
        selection = self._directory / 'subset.txt'
        selection.write_text('${EXPLORATION_QUERY_ROOT}/n.textproto\n$EXPLORATION_QUERY_ROOT/c.textproto\n')
        with mock.patch.dict(os.environ, EXPLORATION_QUERY_ROOT=str(self._directory)):
            queries = self.load(selection)
        self.assertEqual([q.source for q in queries], [str(nitrogen), str(carbon)])

    def test_explicit_name_overrides_file_name_in_legacy_config(self):
        query = self.write_query('n.textproto', 'nitrogen', 'N')
        config = self._directory / 'custom.textproto'
        config.write_text('query { name: "amine_override" fname: "n.textproto" }\n')
        queries = self.load(config)
        self.assertEqual(queries[0].name, 'amine_override')
        self.assertEqual(queries[0].source, str(query))

    def test_duplicate_names_are_rejected(self):
        self.write_query('n.textproto', 'same_name', 'N')
        self.write_query('c.textproto', 'same_name', 'C')
        selection = self._directory / 'queries'
        selection.write_text('n.textproto\nc.textproto\n')
        self.assert_bad(selection, 'Duplicate query names: same_name')

    def test_missing_malformed_and_unnamed_queries(self):
        selection = self._directory / 'queries'
        selection.write_text('missing.textproto\n')
        self.assert_bad(selection, 'Cannot read query file')
        (self._directory / 'missing.textproto').write_text('not a textproto')
        self.assert_bad(selection, 'Cannot parse query file')
        (self._directory / 'missing.textproto').write_text('query { smarts: "N" }')
        self.assert_bad(selection, 'must contain a name')

    def test_compatibility_catalogue_and_subsets(self):
        directory = BIN.parent / 'data/substructure_exploration'
        complete = self.load(directory / 'queries')
        compatibility = self.load(BIN.parent / 'data/substructure_exploration.textproto')
        self.assertEqual([(q.name, q.source) for q in complete],
                         [(q.name, q.source) for q in compatibility])
        functional = self.load(directory / 'functional_groups.txt')
        rings = self.load(directory / 'rings.txt')
        self.assertTrue(functional)
        self.assertTrue(rings)
        self.assertIn('acid', [q.name for q in functional])
        self.assertIn('terminal_aromatic_system', [q.name for q in rings])
        self.assertEqual({q.name for q in functional} | {q.name for q in rings},
                         {q.name for q in complete})

    def test_median_option_analysis_and_html_commands(self):
        self.write_query('n.textproto', 'nitrogen', 'N')
        self.write_query('o.textproto', 'oxygen', 'O')
        selection = self._directory / 'queries'
        selection.write_text('n.textproto\no.textproto\n')
        source = self._directory / 'input.smi'
        source.write_text('CN amine\nCO alcohol\nCC ethane\n')
        activity = self._directory / 'activity.txt'
        activity.write_text('ID activity\namine 3\nalcohol 2\nethane 1\n')
        output = self._directory / 'report.json'
        document = self._directory / 'report.html'
        command = [str(BIN / 'substructure_exploration.sh')]
        env = dict(os.environ, PYTHON=sys.executable)
        analysed = subprocess.run(command + ['analyse', '--config', str(selection),
            '-A', str(activity), '-o', str(output), '--html', str(document),
            '--page_size', '1', str(source)], env=env, capture_output=True, text=True, timeout=20)
        self.assertEqual(analysed.returncode, 0, analysed.stderr)
        self.assertNotIn('median', document.read_text().lower())
        self.assertIn("sortTable('matched_mean')", document.read_text())
        self.assertIn("sortTable('unmatched_mean')", document.read_text())
        self.assertIn('const pageSize = 1;', document.read_text())
        data = json.loads(output.read_text())
        self.assertIn('median', data['global']['summary'])
        self.assertIn('median', data['queries'][0]['matched']['summary'])
        self.assertIn('delta_median', data['queries'][0])
        # Make mean and median ranks deliberately disagree.
        data['queries'][0]['delta_mean'] = 1
        data['queries'][0]['delta_median'] = 100
        data['queries'][1]['delta_mean'] = 100
        data['queries'][1]['delta_median'] = 1
        output.write_text(json.dumps(data))
        original = output.read_bytes()
        for extra, expected in (([], 'oxygen'), (['--sort', 'abs_median'], 'nitrogen')):
            # --median is not given: medians stay hidden even when ranking by them.
            rendered = subprocess.run(command + ['html', str(output), '-o', str(document),
                '--page_size', '0', '--top_n', '1', *extra], env=env,
                capture_output=True, text=True, timeout=20)
            self.assertEqual(rendered.returncode, 0, rendered.stderr)
            text = document.read_text()
            self.assertNotIn('median', text.lower())
            self.assertIn('const pageSize = 0;', text)
            self.assertIn(f'<strong>{expected}</strong>', text)
        self.assertEqual(output.read_bytes(), original)

    def test_selection_cli_from_another_directory(self):
        self.write_query('n.textproto', 'nitrogen', 'N')
        (self._directory / 'queries').write_text('n.textproto\n')
        work = self._directory / 'work'
        work.mkdir()
        (work / 'input.smi').write_text('CN amine\nCC hydrocarbon\n')
        (work / 'activity.txt').write_text('ID activity\namine 3\nhydrocarbon 1\n')
        result = subprocess.run(
            [str(BIN / 'substructure_exploration.sh'), 'analyse', '--config',
             str(self._directory / 'queries'), '-A', 'activity.txt', 'input.smi'],
            cwd=work, env=dict(os.environ, PYTHON=sys.executable),
            capture_output=True, text=True, timeout=20)
        self.assertEqual(result.returncode, 0, result.stderr)
        data = json.loads(result.stdout)
        self.assertEqual(len(data['queries']), 1)
        self.assertEqual(data['queries'][0]['name'], 'nitrogen')
        self.assertEqual(data['queries'][0]['matched']['summary']['n'], 1)


if __name__ == '__main__':
    unittest.main()
