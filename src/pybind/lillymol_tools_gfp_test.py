import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(__file__))

from lillymol import Molecule
from lillymol_tools import GFP, GFPContext, GFPList


def _write_file(contents):
    fd, fname = tempfile.mkstemp(suffix='.gfp')
    with os.fdopen(fd, 'w') as output:
        output.write(contents)
    return fname


def _testdata_file(fname):
    candidates = [
        os.path.join(os.path.dirname(__file__), 'testdata', fname),
    ]
    test_srcdir = os.environ.get('TEST_SRCDIR')
    if test_srcdir:
        candidates.extend([
            os.path.join(test_srcdir, '_main', 'pybind', 'testdata', fname),
            os.path.join(test_srcdir, 'pybind', 'testdata', fname),
        ])
    runfiles_dir = os.environ.get('RUNFILES_DIR')
    if runfiles_dir:
        candidates.extend([
            os.path.join(runfiles_dir, '_main', 'pybind', 'testdata', fname),
            os.path.join(runfiles_dir, 'pybind', 'testdata', fname),
        ])
    for candidate in candidates:
        if os.path.exists(candidate):
            return candidate
    raise FileNotFoundError(fname)

def _mol(smiles):
    mol = Molecule()
    if not mol.build_from_smiles(smiles):
        raise ValueError(smiles)
    return mol


CARBON_GFP = """
$SMI<C>
PCN<Methane>
FCTS<.E..........2;1;1;1;1>
|
$SMI<CC>
PCN<Ethane>
FCTS<.U..........2;1;1;1;1>
|
"""


class TestGFPList(unittest.TestCase):

    def test_read_and_distance(self):
        fname = _write_file(CARBON_GFP)
        try:
            gfp = GFPList.from_file(fname)
            self.assertEqual(len(gfp), 2)
            self.assertEqual(gfp.size(), 2)
            self.assertEqual(gfp.tags(), ['FCTS<'])
            self.assertEqual(gfp.smiles(0), 'C')
            self.assertEqual(gfp.id(1), 'Ethane')
            self.assertAlmostEqual(gfp.distance(0, 0), 0.0, places=6)
            self.assertAlmostEqual(gfp.distance(0, 1), 0.5, places=6)
            self.assertAlmostEqual(gfp.distance(1, 0), 0.5, places=6)
        finally:
            os.remove(fname)

    def test_nearest_neighbours(self):
        fname = _write_file(CARBON_GFP)
        try:
            gfp = GFPList.from_file(fname)
            hits = gfp.nearest_neighbours(0, 1)
            self.assertEqual(len(hits), 1)
            self.assertEqual(hits[0].index, 1)
            self.assertAlmostEqual(hits[0].distance, 0.5, places=6)

            close = gfp.nearest_neighbours_within_distance(0, 0.5)
            self.assertEqual(len(close), 1)
            self.assertEqual(close[0].index, 1)
            self.assertAlmostEqual(close[0].distance, 0.5, places=6)

            none = gfp.nearest_neighbours_within_distance(0, 0.49)
            self.assertEqual(none, [])
        finally:
            os.remove(fname)

    def test_standard_gfp_golden_distances(self):
        gfp = GFPList.from_file(_testdata_file('rand10.standard.gfp'))
        self.assertGreater(len(gfp), 3)
        self.assertEqual(gfp.id(0), 'CHEMBL3460651')
        self.assertEqual(gfp.id(1), 'CHEMBL3460651.a')
        self.assertEqual(gfp.id(3), 'CHEMBL1417367')

        self.assertAlmostEqual(gfp.distance(0, 1), 0.0421, delta=0.0001)
        self.assertAlmostEqual(gfp.distance(1, 0), 0.0421, delta=0.0001)
        self.assertAlmostEqual(gfp.distance(3, 0), 0.499, delta=0.001)
        self.assertAlmostEqual(gfp.distance(0, 3), 0.499, delta=0.001)

    def test_complex_gfp_file(self):
        gfp = GFPList.from_file(_testdata_file('rand10.complex.gfp'))
        self.assertEqual(len(gfp), 10)
        tags = gfp.tags()
        self.assertGreaterEqual(len(tags), 9)
        self.assertIn('FPIW<', tags)
        self.assertIn('FPMK<', tags)
        self.assertIn('MPR<', tags)

        self.assertAlmostEqual(gfp.distance(0, 0), 0.0, places=6)
        d01 = gfp.distance(0, 1)
        self.assertGreaterEqual(d01, 0.0)
        self.assertLessEqual(d01, 1.0)
        self.assertAlmostEqual(d01, gfp.distance(1, 0), places=6)

        hits = gfp.nearest_neighbours(0, 3)
        self.assertEqual(len(hits), 3)
        self.assertNotIn(0, [hit.index for hit in hits])
        self.assertEqual([hit.distance for hit in hits], sorted(hit.distance for hit in hits))

        within = gfp.nearest_neighbours_within_distance(0, 1.0)
        self.assertEqual(len(within), 9)
        self.assertNotIn(0, [hit.index for hit in within])

    def test_use_only_and_use_all(self):
        fname = _write_file(CARBON_GFP)
        try:
            gfp = GFPList.from_file(fname)
            gfp.use_only(['FCTS<'])
            self.assertAlmostEqual(gfp.distance(0, 1), 0.5, places=6)
            gfp.use_all()
            self.assertAlmostEqual(gfp.distance(0, 1), 0.5, places=6)
        finally:
            os.remove(fname)

    def test_generator_specs_match_standard_context(self):
        standard = GFPContext.standard()
        from_specs = GFPContext.from_specs([GFP.iw(), GFP.maccs(), GFP.mpr()])

        self.assertEqual(from_specs.tags(), standard.tags())

        mol = _mol('CCO ethanol')
        fp_standard = standard.fingerprint(mol)
        fp_from_specs = from_specs.fingerprint(mol)
        self.assertAlmostEqual(standard.distance(fp_standard, fp_from_specs),
                               0.0, places=6)

    def test_generator_specs_are_order_independent(self):
        ctx1 = GFPContext.from_specs([GFP.iw(), GFP.maccs(), GFP.mpr()])
        ctx2 = GFPContext.from_specs([GFP.mpr(), GFP.maccs(), GFP.iw()])

        self.assertEqual(ctx1.tags(), ['FPIW<', 'FPMK<', 'FPMK2<', 'MPR<'])
        self.assertEqual(ctx1.tags(), ctx2.tags())

        mol = _mol('CCC propane')
        fp1 = ctx1.fingerprint(mol)
        fp2 = ctx2.fingerprint(mol)
        self.assertAlmostEqual(ctx1.distance(fp1, fp2), 0.0, places=6)

    def test_generator_specs_reject_duplicates(self):
        with self.assertRaises(RuntimeError):
            GFPContext.from_specs([GFP.iw(), GFP.iw()])

    def test_maccs_level2_false(self):
        spec = GFP.maccs(level2=False)
        self.assertEqual(spec.components(), ['FPMK<'])
        self.assertEqual(repr(spec), 'GFP.maccs(level2=False)')

        ctx = GFPContext.from_specs([spec])
        self.assertEqual(ctx.tags(), ['FPMK<'])
        mol = _mol('CCO ethanol')
        fp = ctx.fingerprint(mol)
        self.assertAlmostEqual(ctx.distance(fp, fp), 0.0, places=6)



    def test_ec_generator_spec(self):
        spec = GFP.ec(radius=3, atom_type='UST:AY')
        self.assertEqual(spec.components(), ['NCEC3USTAY<'])
        self.assertEqual(repr(spec), "GFP.ec(radius=3, atom_type='UST:AY')")

        ctx = GFPContext.from_specs([GFP.ec(), spec])
        self.assertEqual(ctx.tags(), ['NCEC3USTAY<', 'NCEC3USTZ<'])
        mol = _mol('CCO ethanol')
        fp = ctx.fingerprint(mol)
        self.assertAlmostEqual(ctx.distance(fp, fp), 0.0, places=6)

    def test_ec_rejects_invalid_specs(self):
        with self.assertRaises(ValueError):
            GFP.ec(radius=-1)
        with self.assertRaises(ValueError):
            GFP.ec(atom_type='')
        with self.assertRaises(RuntimeError):
            GFPContext.from_specs([GFP.ec(radius=3, atom_type='UST:A-Y')])

    def test_alogp_generator_spec(self):
        spec = GFP.alogp(replicates=4)
        self.assertEqual(spec.components(), ['NCALOGP4<'])
        self.assertEqual(repr(spec), 'GFP.alogp(replicates=4)')

        ctx = GFPContext.from_specs([spec])
        self.assertEqual(ctx.tags(), ['NCALOGP4<'])
        mol = _mol('CCO ethanol')
        fp = ctx.fingerprint(mol)
        self.assertAlmostEqual(ctx.distance(fp, fp), 0.0, places=6)

    def test_alogp_rejects_invalid_replicates(self):
        with self.assertRaises(ValueError):
            GFP.alogp(replicates=0)

    def test_standard_context_generates_fingerprints(self):
        ctx = GFPContext.standard()
        self.assertEqual(ctx.tags(), ['FPIW<', 'FPMK<', 'FPMK2<', 'MPR<'])
        self.assertTrue(ctx.can_generate_fingerprints())

        ethane = Molecule()
        self.assertTrue(ethane.build_from_smiles('CC ethane'))
        propane = Molecule()
        self.assertTrue(propane.build_from_smiles('CCC propane'))

        fp1 = ctx.fingerprint(ethane)
        fp2 = ctx.fingerprint(propane)
        self.assertAlmostEqual(ctx.distance(fp1, fp1), 0.0, places=6)
        d12 = ctx.distance(fp1, fp2)
        self.assertGreaterEqual(d12, 0.0)
        self.assertLessEqual(d12, 1.0)
        self.assertAlmostEqual(ctx.distance(fp2, fp1), d12, places=6)

    def test_standard_list_add_and_query_fingerprint(self):
        gfp = GFPList.standard()
        self.assertEqual(gfp.tags(), ['FPIW<', 'FPMK<', 'FPMK2<', 'MPR<'])

        for smiles in ['CC ethane', 'CCC propane', 'CCCC butane']:
            mol = Molecule()
            self.assertTrue(mol.build_from_smiles(smiles))
            gfp.add(mol)

        self.assertEqual(len(gfp), 3)
        self.assertEqual(gfp.id(0), 'ethane')
        self.assertAlmostEqual(gfp.distance(0, 0), 0.0, places=6)

        query = Molecule()
        self.assertTrue(query.build_from_smiles('CCC query'))
        fp = GFPContext.standard().fingerprint(query)
        self.assertAlmostEqual(gfp.distance(fp, 1), 0.0, places=6)

        hits = gfp.nearest_neighbours(fp, 2)
        self.assertEqual(hits[0].index, 1)
        self.assertAlmostEqual(hits[0].distance, 0.0, places=6)


    def test_programmatic_standard_list_matches_stored_file_distances(self):
        stored = GFPList.from_file(_testdata_file('rand10.standard.gfp'))
        molecules = [_mol(stored.smiles(i)) for i in range(len(stored))]
        generated = GFPList.standard_from_molecules(molecules)

        self.assertEqual(len(generated), len(stored))
        self.assertFalse(generated.metadata_stored())
        with self.assertRaises(RuntimeError):
            generated.id(0)
        with self.assertRaises(RuntimeError):
            generated.smiles(0)

        for i in range(len(stored)):
            for j in range(len(stored)):
                self.assertAlmostEqual(generated.distance(i, j),
                                       stored.distance(i, j),
                                       places=6,
                                       msg=f'{i},{j}')

    def test_standard_from_molecules_can_store_metadata(self):
        molecules = [_mol('CC ethane'), _mol('CCC propane'), _mol('CCCC butane')]
        gfp = GFPList.standard_from_molecules(molecules, store_metadata=True)

        self.assertEqual(len(gfp), 3)
        self.assertTrue(gfp.metadata_stored())
        self.assertEqual(gfp.id(0), 'ethane')
        self.assertEqual(gfp.smiles(1), 'CCC')
        self.assertAlmostEqual(gfp.distance(1, 1), 0.0, places=6)

    def test_add_molecules_defaults_to_no_metadata(self):
        molecules = [_mol('CC ethane'), _mol('CCC propane')]
        gfp = GFPList.standard()
        gfp.add_molecules(molecules)

        self.assertEqual(len(gfp), 2)
        self.assertFalse(gfp.metadata_stored())
        self.assertAlmostEqual(gfp.distance(0, 0), 0.0, places=6)
        with self.assertRaises(RuntimeError):
            gfp.id(0)
        with self.assertRaises(RuntimeError):
            gfp.add(_mol('CCCC butane'))

    def test_gfp_failures_raise_exceptions(self):
        mol = Molecule()
        self.assertTrue(mol.build_from_smiles('CC ethane'))

        with self.assertRaises(RuntimeError):
            GFPContext().fingerprint(mol)

        with self.assertRaises(RuntimeError):
            GFPList().add(mol)

        with self.assertRaises(RuntimeError):
            GFPList.from_file('/no/such/gfp/file.gfp')

        gfp = GFPList.standard()
        gfp.add(mol)
        with self.assertRaises(IndexError):
            gfp.distance(0, 1)
        with self.assertRaises(IndexError):
            gfp.smiles(1)
        with self.assertRaises(IndexError):
            gfp.nearest_neighbours(1, 1)
        with self.assertRaises(ValueError):
            gfp.nearest_neighbours_within_distance(0, -1.0)
        with self.assertRaises(RuntimeError):
            gfp.use_only(['NO_SUCH_TAG<'])

        carbon_fname = _write_file(CARBON_GFP)
        try:
            carbon = GFPList.from_file(carbon_fname)
            fp = GFPContext.standard().fingerprint(mol)
            with self.assertRaises(ValueError):
                carbon.distance(fp, 0)
            with self.assertRaises(ValueError):
                carbon.nearest_neighbours(fp, 1)
        finally:
            os.remove(carbon_fname)


    def test_generated_standard_fingerprints_match_stored_file(self):
        stored = GFPList.from_file(_testdata_file('rand10.standard.gfp'))
        self.assertEqual(stored.tags(), ['FPIW<', 'FPMK<', 'FPMK2<', 'MPR<'])

        ctx = GFPContext.standard()
        for i in range(len(stored)):
            mol = Molecule()
            self.assertTrue(mol.build_from_smiles(stored.smiles(i)))
            fp = ctx.fingerprint(mol)
            self.assertAlmostEqual(stored.distance(fp, i), 0.0, places=6,
                                   msg=stored.id(i))

    def test_standard_generation_can_skip_preprocessing(self):
        ctx = GFPContext.standard(preprocess=False)
        self.assertTrue(ctx.can_generate_fingerprints())
        mol = Molecule()
        self.assertTrue(mol.build_from_smiles('CCO ethanol'))
        fp = ctx.fingerprint(mol)
        self.assertAlmostEqual(ctx.distance(fp, fp), 0.0, places=6)

        gfp = GFPList.standard(preprocess=False)
        gfp.add(mol)
        self.assertEqual(len(gfp), 1)
        self.assertAlmostEqual(gfp.distance(0, 0), 0.0, places=6)


if __name__ == '__main__':
    unittest.main()
