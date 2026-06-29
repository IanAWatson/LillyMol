import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(__file__))

from lillymol import *
from lillymol_tools import *


def _skip_guards(test_case):
    """Return True if the test should be skipped due to missing dependencies."""
    try:
        import numpy  # pylint: disable=import-outside-toplevel,unused-import
    except ImportError:
        test_case.skipTest("NumPy is not installed")
    if "LILLYMOL_HOME" not in os.environ:
        test_case.skipTest("LILLYMOL_HOME is not defined")


def _make_molecule(smiles, name):
    mol = Molecule()
    if not mol.build_from_smiles(f"{smiles} {name}"):
        raise ValueError(f"Invalid SMILES: {smiles}")
    return mol


# ---------------------------------------------------------------------------
# Existing single-molecule tests (unchanged)
# ---------------------------------------------------------------------------

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
        self.assertEqual(values[0], 3)   # natoms

        self.assertEqual(names[1], "nrings")
        self.assertEqual(values[1], 0)   # nrings


# ---------------------------------------------------------------------------
# process_list — NumPy array output
# ---------------------------------------------------------------------------

class TestProcessListNumpy(unittest.TestCase):

    def setUp(self):
        _skip_guards(self)
        self._iwdescr = IWDescr()
        self._molecules = [
            _make_molecule("CCO",        "ethanol"),
            _make_molecule("c1ccccc1",   "benzene"),
            _make_molecule("CC(=O)O",    "acetic_acid"),
            _make_molecule("c1ccccc1O",  "phenol"),
        ]

    def test_return_type_is_numpy_array(self):
        import numpy as np
        result = self._iwdescr.process_list(self._molecules)
        self.assertIsInstance(result, np.ndarray)

    def test_shape(self):
        result = self._iwdescr.process_list(self._molecules)
        ndescr = len(self._iwdescr.feature_names())
        self.assertEqual(result.ndim, 2)
        self.assertEqual(result.shape, (len(self._molecules), ndescr))

    def test_dtype_is_float32(self):
        result = self._iwdescr.process_list(self._molecules)
        self.assertEqual(result.dtype.name, "float32")

    def test_row_matches_single_process(self):
        """Each row of process_list must equal the result of process() for that molecule."""
        import numpy as np
        # Re-build molecules so they have not been modified by setUp's process_list call.
        molecules = [
            _make_molecule("CCO",        "ethanol"),
            _make_molecule("c1ccccc1",   "benzene"),
            _make_molecule("CC(=O)O",    "acetic_acid"),
            _make_molecule("c1ccccc1O",  "phenol"),
        ]
        batch = self._iwdescr.process_list(molecules)

        for i, smi_name in enumerate([
            ("CCO",        "ethanol"),
            ("c1ccccc1",   "benzene"),
            ("CC(=O)O",    "acetic_acid"),
            ("c1ccccc1O",  "phenol"),
        ]):
            mol = _make_molecule(*smi_name)
            single = self._iwdescr.process(mol)
            np.testing.assert_array_equal(
                batch[i], single,
                err_msg=f"Row {i} ({smi_name[1]}) differs from process()"
            )

    def test_column_order_matches_feature_names(self):
        """Descriptor values are in the same column order as feature_names()."""
        names = self._iwdescr.feature_names()
        result = self._iwdescr.process_list(self._molecules)
        name_to_col = {n: c for c, n in enumerate(names)}

        # benzene (row 1) has 1 ring
        self.assertEqual(result[1, name_to_col["nrings"]], 1.0)
        # ethanol (row 0) has 0 rings
        self.assertEqual(result[0, name_to_col["nrings"]], 0.0)

    def test_empty_list_returns_zero_rows(self):
        """An empty molecule list should return a (0, ndescr) array, not an error."""
        result = self._iwdescr.process_list([])
        ndescr = len(self._iwdescr.feature_names())
        self.assertEqual(result.shape, (0, ndescr))

    def test_single_molecule_list(self):
        """A one-element list should return a (1, ndescr) array."""
        mol = _make_molecule("CCO", "ethanol")
        result = self._iwdescr.process_list([mol])
        ndescr = len(self._iwdescr.feature_names())
        self.assertEqual(result.shape, (1, ndescr))

    def test_default_does_not_return_dataframe(self):
        """Without as_dataframe=True the result must not be a DataFrame."""
        try:
            import pandas as pd  # pylint: disable=import-outside-toplevel
        except ImportError:
            return  # pandas not installed; nothing to check
        result = self._iwdescr.process_list(self._molecules)
        self.assertNotIsInstance(result, pd.DataFrame)


# ---------------------------------------------------------------------------
# process_list — pandas DataFrame output
# ---------------------------------------------------------------------------

class TestProcessListDataFrame(unittest.TestCase):

    def setUp(self):
        _skip_guards(self)
        try:
            import pandas  # pylint: disable=import-outside-toplevel,unused-import
        except ImportError:
            self.skipTest("pandas is not installed")
        self._iwdescr = IWDescr()
        self._smiles_names = [
            ("CCO",        "ethanol"),
            ("c1ccccc1",   "benzene"),
            ("CC(=O)O",    "acetic_acid"),
            ("c1ccccc1O",  "phenol"),
        ]
        self._molecules = [_make_molecule(s, n) for s, n in self._smiles_names]

    def test_return_type_is_dataframe(self):
        import pandas as pd
        result = self._iwdescr.process_list(self._molecules, as_dataframe=True)
        self.assertIsInstance(result, pd.DataFrame)

    def test_shape(self):
        result = self._iwdescr.process_list(self._molecules, as_dataframe=True)
        ndescr = len(self._iwdescr.feature_names())
        self.assertEqual(result.shape, (len(self._molecules), ndescr))

    def test_column_names_match_feature_names(self):
        names = self._iwdescr.feature_names()
        result = self._iwdescr.process_list(self._molecules, as_dataframe=True)
        self.assertEqual(list(result.columns), names)

    def test_row_index_contains_molecule_names(self):
        result = self._iwdescr.process_list(self._molecules, as_dataframe=True)
        expected_index = [name for _, name in self._smiles_names]
        self.assertEqual(list(result.index), expected_index)

    def test_values_match_numpy_output(self):
        """DataFrame values must be numerically identical to the numpy array output."""
        import numpy as np
        # Rebuild so neither call's molecule-mutation affects the other.
        mols_np = [_make_molecule(s, n) for s, n in self._smiles_names]
        mols_df = [_make_molecule(s, n) for s, n in self._smiles_names]
        arr = self._iwdescr.process_list(mols_np)
        df  = self._iwdescr.process_list(mols_df, as_dataframe=True)
        np.testing.assert_array_equal(arr, df.values)

    def test_lookup_by_column_name(self):
        """DataFrame columns are accessible by descriptor name."""
        result = self._iwdescr.process_list(self._molecules, as_dataframe=True)
        # benzene has 1 ring; ethanol has 0
        self.assertEqual(result.loc["benzene",  "nrings"], 1.0)
        self.assertEqual(result.loc["ethanol",  "nrings"], 0.0)

    def test_lookup_by_row_index(self):
        """DataFrame rows are accessible by molecule name."""
        result = self._iwdescr.process_list(self._molecules, as_dataframe=True)
        row = result.loc["ethanol"]
        self.assertEqual(len(row), len(self._iwdescr.feature_names()))

    def test_dtype_is_float32(self):
        result = self._iwdescr.process_list(self._molecules, as_dataframe=True)
        self.assertEqual(result["natoms"].dtype.name, "float32")

    def test_empty_list_returns_empty_dataframe(self):
        """An empty molecule list should return a DataFrame with 0 rows and correct columns."""
        result = self._iwdescr.process_list([], as_dataframe=True)
        names = self._iwdescr.feature_names()
        self.assertEqual(result.shape[0], 0)
        self.assertEqual(list(result.columns), names)

    def test_single_molecule_list(self):
        """A one-element list should return a single-row DataFrame."""
        mol = _make_molecule("CCO", "ethanol")
        result = self._iwdescr.process_list([mol], as_dataframe=True)
        self.assertEqual(result.shape[0], 1)
        self.assertEqual(result.index[0], "ethanol")

    def test_false_flag_returns_numpy_not_dataframe(self):
        """Explicitly passing as_dataframe=False must return a numpy array."""
        import numpy as np
        import pandas as pd
        result = self._iwdescr.process_list(self._molecules, as_dataframe=False)
        self.assertIsInstance(result, np.ndarray)
        self.assertNotIsInstance(result, pd.DataFrame)


if __name__ == "__main__":
    unittest.main()
