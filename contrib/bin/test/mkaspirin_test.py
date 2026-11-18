#!/usr/bin/env python3
"""Run through LillyMol's run_python.sh with staged bindings; RDKit is optional."""
from pathlib import Path
import subprocess
import sys
import unittest

BIN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BIN))
import mkaspirin
from lillymol import MolFromSmiles, SubstructureQuery

PANEL = [
    (mkaspirin.ASPIRIN, True),
    ('O=C(O)c1ccccc1O', False),
    ('CC(=O)Oc1ccccc1C(=O)OC', False),
    ('CC(=O)Oc1ccccc1', False),
    ('CC(=O)Oc1cccc(C(=O)O)c1', False),
    ('CC(=O)Oc1ccc(C(=O)O)cc1', False),
    ('CCC(=O)Oc1ccccc1C(=O)O', False),
    ('CC(=O)Oc1ccccc1C(=O)[O-]', False),
    (mkaspirin.ASPIRIN + '.[Na+]', True),
    ('[13CH3]C(=O)Oc1ccccc1C(=O)O', False),
    ('CC(=O)Oc1ccc(C)cc1C(=O)O', False),
]


class GeneratorTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = mkaspirin.AspirinModel()
        cls.molecules = [MolFromSmiles(smiles) for smiles, _ in PANEL]
        cls.reordered = [MolFromSmiles(cls.model.molecule.random_smiles()) for _ in range(5)]

    def generator(self, **kwargs):
        return mkaspirin.Generator(self.model, **kwargs)

    def query(self, smarts):
        query = SubstructureQuery()
        self.assertTrue(query.build_from_smarts(smarts))
        return query

    def check_matches(self, generator):
        smarts = generator.smarts()
        self.assertEqual(len(smarts), generator.size)
        query = self.query(smarts)
        for mol, (source, expected) in zip(self.molecules, PANEL):
            self.assertEqual(query.substructure_search(mol) > 0, expected, source)
        for mol in self.reordered:
            self.assertGreater(query.substructure_search(mol), 0)

    def test_every_context_is_rooted_at_the_correct_source_atom(self):
        generator = self.generator()
        for atom in range(13):
            for context in mkaspirin.CONTEXTS:
                with self.subTest(atom=atom, context=context):
                    generator.parts, _, _ = generator._instantiate(atom, context, 0)
                    query = self.query(generator.smarts())
                    results = query.substructure_search_embeddings(self.model.molecule)
                    self.assertIn(atom, {embedding[0] for embedding in results})

    def test_baseline_and_expansions_preserve_matches(self):
        for seed in (1, 2026):
            for extensions in (False, True):
                for alternatives in (1, 2):
                    with self.subTest(seed=seed, extensions=extensions, alternatives=alternatives):
                        generator = self.generator(seed=seed, alternatives=alternatives,
                                                   lillymol_extensions=extensions)
                        self.check_matches(generator)
                        for step in range(12):
                            generator.expand()
                            if step in (0, 5, 11):
                                self.check_matches(generator)
                        text = generator.smarts()
                        self.assertEqual('/IWfss' in text, extensions)
                        self.assertEqual('>' in text or '<' in text, extensions)

    def test_seed_is_reproducible(self):
        a, b, c = [self.generator(seed=seed) for seed in (6, 6, 7)]
        for _ in range(25):
            for generator in (a, b, c):
                generator.expand()
        self.assertEqual(a.smarts(), b.smarts())
        self.assertNotEqual(a.smarts(), c.smarts())

    def test_depth_limit_exhaustion(self):
        generator = self.generator(max_depth=1, contexts=('molecule',))
        for _ in range(13):
            generator.expand()
        self.assertEqual(generator.depth, 1)
        self.assertEqual(len(generator.smarts()), generator.size)
        with self.assertRaisesRegex(ValueError, 'No expandable atoms'):
            generator.expand()

    def test_serialization_does_not_depend_on_python_recursion_limit(self):
        generator = self.generator(max_depth=1200, contexts=('atom',))
        # Restrict the initial expansion sites to force a single deep chain.
        generator.expandable = generator.expandable[:1]
        for _ in range(1100):
            generator.expand()
        self.assertEqual(generator.depth, 1100)
        self.assertEqual(len(generator.smarts()), generator.size)

    def test_portable_patterns_also_match_in_rdkit(self):
        try:
            from rdkit import Chem
        except ImportError:
            self.skipTest('RDKit is not installed')
        molecules = [Chem.MolFromSmiles(smiles) for smiles, _ in PANEL]
        for seed in (3, 42):
            generator = self.generator(seed=seed, alternatives=2)
            for step in range(17):
                if step in (0, 4, 16):
                    query = Chem.MolFromSmarts(generator.smarts())
                    self.assertIsNotNone(query)
                    for molecule, (source, expected) in zip(molecules, PANEL):
                        self.assertEqual(molecule.HasSubstructMatch(query), expected, source)
                generator.expand()


class CommandLineTest(unittest.TestCase):
    def run_script(self, *args):
        return subprocess.run([sys.executable, str(BIN/'mkaspirin.py'), *args],
                              capture_output=True, text=True, timeout=30)

    def test_output_counts_and_last_only(self):
        all_rows = self.run_script('--nsmarts', '5', '--seed', '19')
        last = self.run_script('--nsmarts', '5', '--seed', '19', '--write_last')
        self.assertEqual(all_rows.returncode, 0, all_rows.stderr)
        self.assertEqual(last.returncode, 0, last.stderr)
        rows = all_rows.stdout.splitlines()
        self.assertEqual(len(rows), 5)
        self.assertEqual(last.stdout.strip(), rows[-1])
        self.assertEqual([row.split()[1] for row in rows], ['#1', '#2', '#3', '#4', '#5'])

    def test_target_size_and_extension_flag(self):
        result = self.run_script('--target_bytes', '20000', '--write_last', '--lillymol_extensions')
        self.assertEqual(result.returncode, 0, result.stderr)
        smarts = result.stdout.split()[0]
        self.assertGreaterEqual(len(smarts), 20000)
        self.assertIn('/IWfss', smarts)
        self.assertIn(f'{len(smarts)} SMARTS bytes', result.stderr)
        baseline = self.run_script('--nsmarts', '0', '--no_lillymol_extensions')
        self.assertEqual(baseline.returncode, 0, baseline.stderr)
        self.assertEqual(len(baseline.stdout.splitlines()), 1)
        self.assertNotIn('/IW', baseline.stdout)

    def test_bad_arguments_and_insufficient_depth_fail(self):
        for args in (['--nsmarts', '-1'], ['--target_bytes', '0'],
                     ['--alternatives', '0'], ['--contexts', 'unknown'],
                     ['--nsmarts', '1', '--target_bytes', '1000'],
                     ['--lillymol_extensions', '--no_lillymol_extensions'],
                     ['--max_depth', '0', '--nsmarts', '1', '--write_last']):
            with self.subTest(args=args):
                result = self.run_script(*args)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(result.stdout, '')
                self.assertNotIn('Traceback', result.stderr)


if __name__ == '__main__':
    unittest.main()
