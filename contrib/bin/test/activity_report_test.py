"""Standard-library-only tests for shared ranking and fragment HTML output."""
import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

BIN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BIN))
from activity_report import build_html, ranked_rows


def row(name, effect, support=2, absent=2):
    return dict(name=name, source=name, delta_median=effect, delta_mean=effect,
                matched=dict(summary=dict(n=support, median=effect, mean=effect), histogram=[0, support]),
                unmatched=dict(summary=dict(n=absent, median=0, mean=0), histogram=[absent, 0]),
                examples=[dict(identifier='<parent>', smiles='C&N', activity=1)])


def report():
    return dict(metadata=dict(feature_kind='fragment', queries=4, activity='activity.txt',
                              molecules_read_with_activity=4),
                **{'global': dict(summary=dict(median=0), histogram=[2, 2])},
                queries=[row('positive', 4), row('negative', -8),
                         row('sparse', 10, support=1), row('universal', 99, absent=0)])


class ActivityReportTest(unittest.TestCase):
    def test_ranking_support_and_signed_directions(self):
        data = report()
        self.assertEqual([r['name'] for r in ranked_rows(data, top_n=2)], ['sparse', 'negative'])
        self.assertEqual([r['name'] for r in ranked_rows(data, min_matches=2)], ['negative', 'positive'])
        self.assertEqual([r['name'] for r in ranked_rows(data, sort_by='median', top_n=1)], ['sparse'])
        self.assertEqual([r['name'] for r in ranked_rows(data, sort_by='median', ascending=True, top_n=1)], ['negative'])
        self.assertEqual(len(ranked_rows(data, top_n=0, min_unmatched=0)), 4)
        self.assertEqual(len(data['queries']), 4)

    def test_ties_and_missing_deltas(self):
        data = report()
        data['queries'] = [row('Z', -5), row('A', 5), row('missing', 0)]
        data['queries'][-1]['delta_median'] = None
        self.assertEqual([r['name'] for r in ranked_rows(data)], ['A', 'Z', 'missing'])
        with self.assertRaises(ValueError):
            ranked_rows(data, top_n=-1)

    def test_fragment_and_dictionary_rendering(self):
        data = report()
        html = build_html(data, '<title>', rows=ranked_rows(data, top_n=1))
        self.assertIn('1 of 4 fragments', html)
        self.assertIn('<th>Fragment</th>', html)
        self.assertNotIn('>Multiple</th>', html)
        self.assertIn('&lt;parent&gt;', html)
        self.assertIn('C&amp;N', html)
        self.assertIn('&lt;title&gt;', html)
        dictionary = copy.deepcopy(data)
        dictionary['metadata'].pop('feature_kind')
        for item in dictionary['queries']:
            item['multiplicity'] = dict(molecules_with_multiple_matches=0, max_matches=1)
        html = build_html(dictionary, 'Queries')
        self.assertIn('<th>Query</th>', html)
        self.assertIn('>Multiple</th>', html)

    def test_cli_with_no_site_packages_stdin_and_top_n(self):
        with tempfile.TemporaryDirectory(prefix='activity report ') as tmp:
            root = Path(tmp)
            data = report()
            data['queries'] = [row(f'fragment{i}', i) for i in range(101)]
            data['metadata']['queries'] = 101
            source = root / 'report.json'
            source.write_text(json.dumps(data))
            command = [sys.executable, '-S', str(BIN / 'dicer_fragment_activity_html.py')]
            result = subprocess.run(command + [str(source), '-o', '-'], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('100 of 101 fragments', result.stdout)
            self.assertNotIn('<strong>fragment0</strong>', result.stdout)
            result = subprocess.run(command + ['-', '-o', '-', '--top_n', '0'],
                                    input=source.read_text(), capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('101 of 101 fragments', result.stdout)
            before = source.read_bytes()
            result = subprocess.run(command + [str(source), '-o', str(source)], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(source.read_bytes(), before)
            alias = root / 'hardlink.json'
            alias.hardlink_to(source)
            result = subprocess.run(command + [str(source), '-o', str(alias)], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(source.read_bytes(), before)


if __name__ == '__main__':
    unittest.main()
