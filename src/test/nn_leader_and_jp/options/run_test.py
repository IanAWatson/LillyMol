#!/usr/bin/env python3
"""Regression checks for supported clustering options."""
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

EXE = str(Path(sys.argv[1]).resolve())


class Options(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.input = Path(self.tmp.name) / 'input.nn'
        records = []
        distances = [[0, 0.1, 0.8], [0.1, 0, 0.8], [0.8, 0.8, 0]]
        for i in range(3):
            lines = ['$SMI<C>', f'PCN<m{i} 2>']
            for j in sorted(range(3), key=lambda j: distances[i][j]):
                if i != j:
                    lines.extend([f'PCN<m{j}>', f'DIST<{distances[i][j]}>'])
            records.append('\n'.join([*lines, '|', '']))
        self.input.write_text(''.join(records))

    def run_tool(self, *args, success=True):
        result = subprocess.run([EXE, *map(str, args), str(self.input)],
                                capture_output=True, text=True, timeout=15)
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertGreater(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout, '')
        return result

    def test_removed_threshold_option(self):
        for mode in ([], ['-T']):
            result = self.run_tool(*mode, '-H', 'THR<', success=False)
            self.assertIn('Unrecognised', result.stderr)

    def test_removed_random_leader_option(self):
        for mode in ([], ['-T']):
            result = self.run_tool(*mode, '-t', '0.2', '-e', success=False)
            self.assertIn('Unrecognised', result.stderr)

    def test_removed_score_option(self):
        for spec in ('SCORE<', 'col=2', 'FILE=missing.scores'):
            with self.subTest(spec=spec):
                result = self.run_tool('-t', '0.9', '-S', spec, success=False)
                self.assertIn('Unrecognised', result.stderr)

    def priority_records(self, order):
        # m3 is isolated; the others form a clique. Names and SCORE tags
        # deliberately put the largest values on later clique members.
        scores = {'m3': 1000, 'm2': 1, 'm0': 900, 'm1': 500}
        records = []
        for name in order:
            lines = ['$SMI<C>', f'PCN<{name} {scores[name]}>',
                     f'SCORE<{scores[name]}>']
            if name != 'm3':
                for nbr in ('m0', 'm1', 'm2'):
                    if nbr != name:
                        lines.extend([f'PCN<{nbr}>', 'DIST<0.1>'])
            records.append('\n'.join([*lines, '|', '']))
        return ''.join(records)

    def leaders(self, output):
        return [re.search(r'^PCN<(m\d+)', c, re.MULTILINE).group(1)
                for c in output.split('|') if c.strip()]

    def test_input_priority(self):
        for first, last in (('m2', 'm1'), ('m1', 'm2')):
            self.input.write_text(self.priority_records(['m3', first, 'm0', last]))
            # Leader respects input priority even when its first choice is
            # isolated. Butina first maximises degree, then uses input order.
            self.assertEqual(self.leaders(self.run_tool('-h', '-t', '0.2').stdout),
                             ['m3', first])
            self.assertEqual(self.leaders(self.run_tool('-T', '-h', '-t', '0.2').stdout),
                             [first, 'm3'])

    def test_input_priority_across_files(self):
        first = Path(self.tmp.name) / 'first.nn'
        first.write_text(self.priority_records(['m3', 'm2']))
        self.input.write_text(self.priority_records(['m0', 'm1']))
        self.assertEqual(self.leaders(self.run_tool('-s', '4', '-h', '-t', '0.2', first).stdout),
                         ['m3', 'm2'])
        self.assertEqual(self.leaders(self.run_tool('-s', '4', '-T', '-h', '-t', '0.2', first).stdout),
                         ['m2', 'm3'])

    def complete_graph(self, names, tag='PCN<'):
        records = []
        for name in names:
            lines = ['$SMI<C>', f'{tag}{name}>']
            for neighbour in names:
                if neighbour != name:
                    lines.extend([f'{tag}{neighbour.split()[0]}>', 'DIST<0.1>'])
            records.append('\n'.join([*lines, '|', '']))
        self.input.write_text(''.join(records))

    def cluster_ids(self, output, tag='PCN<'):
        return [re.findall(r'^' + re.escape(tag) + r'([^>]+)>', c, re.MULTILINE)
                for c in output.split('|') if c.strip()]

    def test_single_character_identifiers(self):
        for tag in ('PCN<', 'ID<', 'MOLECULE_ID<'):
            for mode in (['-h', '-t', '0.2'], ['-T', '-h', '-t', '0.2'], ['-J', '.,1']):
                with self.subTest(tag=tag, mode=mode):
                    self.complete_graph(['a', 'b', 'c'], tag)
                    result = self.run_tool(*mode, '-I', tag)
                    self.assertEqual(self.cluster_ids(result.stdout, tag), [['a', 'b', 'c']])

    def test_limit_from_entire_first_token(self):
        for tag in ('PCN<', 'ID<', 'MOLECULE_ID<'):
            for mode in ([], ['-T']):
                for names, expected in ((['12', '13', '14', '15'], [['12', '13', '14', '15']]),
                                        (['1', '2', '3', '4'], [['1'], ['2', '3'], ['4']])):
                    with self.subTest(tag=tag, mode=mode, names=names):
                        self.complete_graph(names, tag)
                        result = self.run_tool(*mode, '-h', '-t', '0.2', '-I', tag,
                                               '-M', 'col=1')
                        self.assertEqual(self.cluster_ids(result.stdout, tag), expected)

    def test_tag_cluster_limit_rejected(self):
        result = self.run_tool('-t', '0.9', '-M', 'MAX<', success=False)
        self.assertIn('-M', result.stderr)
        self.assertIn('col=<col>', result.stderr)

    def test_invalid_limit_columns_rejected(self):
        for column in ('', '0', '-1', 'text'):
            with self.subTest(column=column):
                self.run_tool('-t', '0.9', '-M', f'col={column}', success=False)

    def test_column_cluster_limit(self):
        for mode in ([], ['-T']):
            result = self.run_tool(*mode, '-h', '-t', '0.9', '-M', 'col=2')
            clusters = [c for c in result.stdout.split('|') if c.strip()]
            members = [re.findall(r'^PCN<(m\d+)', c, re.MULTILINE) for c in clusters]
            self.assertEqual(sorted(map(len, members)), [1, 2])
            self.assertEqual(sorted(m for cluster in members for m in cluster),
                             ['m0', 'm1', 'm2'])

    def test_threshold_required_for_leader_and_butina(self):
        for mode in ([], ['-T']):
            result = self.run_tool(*mode, success=False)
            self.assertIn('Threshold distance must be specified via -t', result.stderr)

    def endpoint_graph(self):
        # The zero-distance pair is joined at 0; atom m2 is joined only at 1.
        # m3 has no listed neighbours, so it remains isolated even at 1.
        edges = {'m0': [('m1', '0'), ('m2', '1')],
                 'm1': [('m0', '0'), ('m2', '1')],
                 'm2': [('m0', '1'), ('m1', '1')], 'm3': []}
        records = []
        for name, neighbours in edges.items():
            lines = ['$SMI<C>', f'PCN<{name}>']
            for neighbour, distance in neighbours:
                lines.extend([f'PCN<{neighbour}>', f'DIST<{distance}>'])
            records.append('\n'.join([*lines, '|', '']))
        self.input.write_text(''.join(records))

    def test_inclusive_threshold_endpoints(self):
        self.endpoint_graph()
        for mode in ([], ['-T']):
            expected = {}
            for threshold, clusters in (('0', [['m0', 'm1'], ['m2'], ['m3']]),
                                        ('1', [['m0', 'm1', 'm2'], ['m3']])):
                result = self.run_tool(*mode, '-h', '-t', threshold)
                self.assertEqual(self.cluster_ids(result.stdout), clusters)
                expected[threshold] = result.stdout
            for args in (['-t', '0,1,0'], ['-t', '0', '-t', '1', '-t', '0']):
                stem = Path(self.tmp.name) / 'endpoints'
                self.run_tool(*mode, '-h', *args, '-B', stem)
                for i, threshold in enumerate(('0', '1', '0')):
                    self.assertEqual(Path(f'{stem}{i}.ldr').read_text(), expected[threshold])

    def test_invalid_thresholds_rejected_in_all_forms(self):
        for mode in ([], ['-T']):
            for value in ('-0.0001', '1.0001', 'text', 'nan', 'inf'):
                for args in (['-t', value], ['-t', f'0,{value},1'],
                             ['-t', '0', '-t', value, '-t', '1']):
                    with self.subTest(mode=mode, args=args):
                        result = self.run_tool(*mode, *args, success=False)
                        self.assertIn('[0, 1]', result.stderr)

    def test_empty_threshold_value_rejected(self):
        for mode in ([], ['-T']):
            for args in (['-t', ''], ['-t', '0', '-t', '', '-t', '1']):
                result = self.run_tool(*mode, *args, success=False)
                self.assertIn('[0, 1]', result.stderr)

    def test_multiple_global_thresholds(self):
        for mode in ([], ['-T']):
            for thresholds in (['-t', '0.2,0.9'], ['-t', '0.2', '-t', '0.9']):
                stem = Path(self.tmp.name) / 'clusters'
                self.run_tool(*mode, '-h', *thresholds, '-B', stem)
                for i, threshold in enumerate(('0.2', '0.9')):
                    expected = self.run_tool(*mode, '-h', '-t', threshold).stdout
                    self.assertEqual(Path(f'{stem}{i}.ldr').read_text(), expected)


if __name__ == '__main__':
    result = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Options))
    if not result.wasSuccessful():
        sys.exit(1)
    print('nn_leader_and_jp option checks passed')
