"""Run with PYTHON=/path/to/matching/python $LILLYMOL_HOME/run_python.sh this_file."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

BIN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BIN))
import substructure_exploration as exploration
from lillymol import MolFromSmiles, Standardise


class CatalogueTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.queries = exploration.build_queries(exploration.read_config(
            str(BIN.parent / 'data/substructure_exploration.textproto')))
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

    def test_fluorine_partition(self):
        for smiles, cf3, other in [
            ('FC(F)(F)c1ccc(F)cc1', 1, 1), ('CC(F)F', 0, 2),
            ('FC(F)(F)F', 0, 4), ('CC(F)(F)F', 1, 0),
        ]:
            with self.subTest(smiles=smiles):
                self.assertEqual(self.count('cf3', smiles), cf3)
                self.assertEqual(self.count('fluorine_not_cf3', smiles), other)
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
            self.assertIn('<html', (work / 'report.html').read_text())
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



if __name__ == '__main__':
    unittest.main()
