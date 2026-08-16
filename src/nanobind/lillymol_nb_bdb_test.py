import os
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(__file__))


def _runfile_path(path):
    runfiles_dir = os.environ.get("RUNFILES_DIR")
    if runfiles_dir:
        candidate = os.path.join(runfiles_dir, "_main", path)
        if os.path.exists(candidate):
            return candidate

    candidate = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", path))
    if os.path.exists(candidate):
        return candidate

    candidate = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "bazel-bin", path))
    if os.path.exists(candidate):
        return candidate

    return None

import lillymol_nb
import lillymol_nb_bdb


class TestNanobindBdb(unittest.TestCase):

    def test_import_and_construct(self):
        self.assertIsNotNone(lillymol_nb_bdb.Selimsteg())
        precedent = lillymol_nb_bdb.SyntheticPrecedentDatabases()
        self.assertIn("Synthetic Precedent database", repr(precedent))
        self.assertTrue(precedent.set_max_radius(2))

    def test_open_missing_database_fails(self):
        lookup = lillymol_nb_bdb.Selimsteg()
        self.assertFalse(lookup.open_database("/no/such/berkeleydb/file.bdb"))

    def test_selimsteg_lookup(self):
        loader = _runfile_path("BerkeleyDB/iwbdb_load")
        if loader is None:
            self.skipTest("iwbdb_load not available in runfiles")

        with tempfile.TemporaryDirectory() as tmpdir:
            input_fname = os.path.join(tmpdir, "id_smiles.txt")
            dbname = os.path.join(tmpdir, "id_smiles.bdb")
            with open(input_fname, "w", encoding="utf-8") as writer:
                writer.write("mol1 CCO\n")
                writer.write("mol2 c1ccccc1\n")

            subprocess.run([loader, "-d", dbname, input_fname], check=True)

            lookup = lillymol_nb_bdb.Selimsteg()
            self.assertTrue(lookup.open_database(dbname))
            self.assertEqual(lookup.get_smiles("mol1"), "CCO")
            self.assertIsNone(lookup.get_smiles("missing"))

            mol = lookup.get_molecule("mol2")
            self.assertIsNotNone(mol)
            self.assertEqual(mol.unique_smiles(), "c1ccccc1")
            self.assertEqual(mol.name(), "mol2")

            mols = lookup.get_molecules(["mol1", "missing", "mol2"])
            self.assertEqual(len(mols), 3)
            self.assertEqual(mols[0].unique_smiles(), "OCC")
            self.assertEqual(mols[1].natoms(), 0)
            self.assertEqual(mols[2].unique_smiles(), "c1ccccc1")



if __name__ == "__main__":
    unittest.main()
