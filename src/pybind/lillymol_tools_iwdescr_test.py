import os
import unittest

# Codex made these the import path,
# from pybind import lillymol
# from pybind import lillymol_tools
# but then a molecule is m = lillymol.Molecule
# and iwdescr = lillymol_tools.IWDescr
# And it was able to make that work from a bazel test.
# TODO:ianwatson investigate sometime
from lillymol import *
from lillymol_tools import *


class TestIWDescr(unittest.TestCase):

  def test_all_descriptors(self):
    try:
      import numpy  # pylint: disable=import-outside-toplevel,unused-import
    except ImportError:
      self.skipTest("NumPy is not installed")

    if "LILLYMOL_HOME" not in os.environ:
      self.skipTest("LILLYMOL_HOME is not defined")

    mol = Molecule()
    self.assertTrue(mol.build_from_smiles("CCO"))

    iwdescr = IWDescr()
    names = iwdescr.feature_names()
    values = iwdescr.process(mol)

    self.assertGreater(len(names), 0)
    self.assertEqual(len(names), len(set(names)))
    self.assertTrue(all(names))
    self.assertEqual(values.ndim, 1)
    self.assertEqual(values.dtype.name, "float32")
    self.assertEqual(values.size, len(names))

    self.assertEqual(names[0], "natoms")
    self.assertEqual(values[0], 3)  # natoms
    self.assertEqual(names[1], "nrings")
    self.assertEqual(values[1], 0)  # nrings


if __name__ == "__main__":
  unittest.main()
