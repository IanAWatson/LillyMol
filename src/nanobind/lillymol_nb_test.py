import atexit
import copy
import math
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(__file__))

import lillymol_nb


def _trace_process_exit():
    if os.environ.get("LILLYMOL_NB_TEST_TRACE"):
        print("ATEXIT unittest complete", file=sys.stderr, flush=True)


atexit.register(_trace_process_exit)


class LillyMolNanobindTestCase(unittest.TestCase):

    def setUp(self):
        if os.environ.get("LILLYMOL_NB_TEST_TRACE"):
            print(f"START {self.id()}", file=sys.stderr, flush=True)

    def tearDown(self):
        if os.environ.get("LILLYMOL_NB_TEST_TRACE"):
            print(f"END {self.id()}", file=sys.stderr, flush=True)

SMILES = """C(=O)NCCCC CHEMBL45466
C(O)(=O)CCC(O)=O CHEMBL1200345
N1=NC(=CS1)C(O)=O CHEMBL247337
NOC1CCC(N)CC1 CHEMBL1213430
N1(=C(C)CCC1(C)C)=O CHEMBL325242
O=C(NC)[C@H]1NC(=O)CC1 CHEMBL1892080
C1=CC=C2C(=CC=N2)N1C CHEMBL593929
C1=CC=C2C(=C1)NC(N)S2 CHEMBL568765
C1C(N)CC1(O)P(O)(C)=O CHEMBL509338
N1(C(C#N)C1)C(=O)NCC CHEMBL150159
"""


class TestNanobindMolecule(LillyMolNanobindTestCase):

    def test_build_from_smiles(self):
        mol = lillymol_nb.Molecule()
        self.assertTrue(mol.build_from_smiles("CCO ethanol"))
        self.assertEqual(mol.natoms(), 3)
        self.assertEqual(mol.nedges(), 2)
        self.assertEqual(mol.name(), "ethanol")
        self.assertEqual(mol.smiles(), "CCO")
        self.assertEqual(mol.unique_smiles(), "OCC")
        self.assertEqual(mol.molecular_formula(), "C2OH6")

    def test_mol_from_smiles(self):
        mol = lillymol_nb.MolFromSmiles("c1ccccc1 benzene")
        self.assertIsNotNone(mol)
        self.assertEqual(mol.natoms(), 6)
        self.assertEqual(mol.nrings(), 1)
        self.assertEqual(mol.name(), "benzene")
        self.assertIsNone(lillymol_nb.MolFromSmiles("["))

    def test_lillymol_from_smiles_and_batch_smiles(self):
        mol = lillymol_nb.LillyMolFromSmiles("CCO ethanol")
        self.assertIsNotNone(mol)
        self.assertEqual(mol.name(), "ethanol")

        molecules = lillymol_nb.MolFromSmiles(["C methane", "CC ethane", "["])
        self.assertEqual(len(molecules), 3)
        self.assertEqual(molecules[0].name(), "methane")
        self.assertEqual(molecules[1].natoms(), 2)
        self.assertEqual(molecules[2].natoms(), 0)

    def test_molecule_equality(self):
        mol1 = lillymol_nb.MolFromSmiles("CCO ethanol")
        mol2 = lillymol_nb.MolFromSmiles("CCO other")
        mol3 = lillymol_nb.MolFromSmiles("CCN ethylamine")
        self.assertEqual(mol1, mol2)
        self.assertNotEqual(mol1, mol3)
        mol2.set_atomic_number(2, 7)
        self.assertEqual(mol2, mol3)

    def test_reader_open_and_next(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, "input.smi")
            with open(fname, "w") as writer:
                writer.write(SMILES)

            reader = lillymol_nb.Reader()
            self.assertTrue(reader.open(fname))
            mol = reader.next()
            self.assertIsNotNone(mol)
            self.assertEqual(mol.name(), "CHEMBL45466")
            self.assertEqual(reader.molecules_read(), 1)

    def test_reader_bad_suffix_and_explicit_type(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, "input.foo")
            with open(fname, "w") as writer:
                writer.write(SMILES)

            reader = lillymol_nb.Reader()
            self.assertFalse(reader.open(fname))
            self.assertTrue(reader.open(fname, lillymol_nb.FileType.SMI))
            self.assertEqual(reader.next().name(), "CHEMBL45466")

    def test_reader_iter(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, "input.smi")
            with open(fname, "w") as writer:
                writer.write(SMILES)

            reader = lillymol_nb.Reader()
            self.assertTrue(reader.open(fname, lillymol_nb.FileType.SMI))
            names = [mol.name() for mol in reader]
            self.assertEqual(len(names), 10)
            self.assertEqual(reader.molecules_read(), 10)

    def test_slurp(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, "input.smi")
            with open(fname, "w") as writer:
                writer.write(SMILES)

            molecules = lillymol_nb.slurp(fname)
            self.assertIsNotNone(molecules)
            self.assertEqual(len(molecules), 10)
            self.assertEqual(molecules[0].name(), "CHEMBL45466")
            self.assertEqual([mol.name() for mol in molecules[:2]], ["CHEMBL45466", "CHEMBL1200345"])

    def test_global_sdf_option_helpers_are_available(self):
        self.assertTrue(lillymol_nb.set_sdf_identifier("idnumber"))
        lillymol_nb.set_prepend_sdfid(True)
        lillymol_nb.set_allsdfid(False)
        lillymol_nb.set_sdf_tags_to_json(False)
        lillymol_nb.set_firstsdftag(False)
        lillymol_nb.set_ignore_bad_m(False)
        lillymol_nb.set_mdlquiet(False)
        lillymol_nb.set_allow_deuterium(False)
        lillymol_nb.set_allow_tritium(False)

    def test_context_aliases(self):
        self.assertIs(lillymol_nb.ReaderContext, lillymol_nb.MolReaderContext)
        self.assertIs(lillymol_nb.ContextWriter, lillymol_nb.MolWriterContext)

    def test_mol_reader_context(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, "input.smi")
            with open(fname, "w") as writer:
                writer.write(SMILES)

            with lillymol_nb.MolReaderContext(fname, lillymol_nb.FileType.SMI) as reader:
                names = [mol.name() for mol in reader]
            self.assertEqual(len(names), 10)
            self.assertEqual(names[0], "CHEMBL45466")

    def test_mol_reader_context_preprocessing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, "input.smi")
            with open(fname, "w") as writer:
                writer.write("CC.O mixture\n")

            reader = lillymol_nb.MolReaderContext(fname, largest_fragment=True)
            mol = reader.next()
            self.assertIsNotNone(mol)
            self.assertEqual(mol.natoms(), 2)
            self.assertIsNone(reader.next())

    def test_retained_sdf_text_info_and_tags(self):
        sdf = """ethanol
  LillyMol          2D

  3  2  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
>  <ID Number>
CHEMBL1

>  <Comment>
first line
second line

$$$$
"""
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, "input.sdf")
            with open(fname, "w") as writer:
                writer.write(sdf)

            reader = lillymol_nb.MolReaderContext(
                fname, lillymol_nb.FileType.SDF, keep_sdf_tags=True)
            mol = reader.next()
            self.assertIsNotNone(mol)
            self.assertGreater(mol.number_records_text_info(), 0)
            self.assertIn(">  <ID Number>", mol.text_info())
            self.assertEqual(mol.sdf_tags()["ID_Number"], "CHEMBL1")
            self.assertEqual(mol.sdf_tags()["Comment"], "first line\nsecond line")

    def test_writer_and_context_writer(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            stem = os.path.join(tmpdir, "written")
            mol = lillymol_nb.MolFromSmiles("CCO ethanol")
            writer = lillymol_nb.Writer()
            self.assertTrue(writer.add_output_type(lillymol_nb.FileType.SMI))
            self.assertTrue(writer.new_stem(stem))
            self.assertTrue(writer.write(mol))
            writer.close()

            reader = lillymol_nb.Reader()
            self.assertTrue(reader.open(stem + ".smi"))
            self.assertEqual(reader.next().name(), "ethanol")

            stem = os.path.join(tmpdir, "context")
            with lillymol_nb.MolWriterContext(stem, lillymol_nb.FileType.SMI) as writer:
                self.assertTrue(writer.write(mol))

            reader = lillymol_nb.Reader()
            self.assertTrue(reader.open(stem + ".smi"))
            self.assertEqual(reader.next().smiles(), "CCO")

    def test_set_name(self):
        mol = lillymol_nb.Molecule()
        self.assertTrue(mol.build_from_smiles("C methane"))
        mol.set_name("renamed")
        self.assertEqual(mol.name(), "renamed")

    def test_rdkit_style_molecule_access_aliases(self):
        mol = lillymol_nb.MolFromSmiles("[H]OC ethanol")
        self.assertEqual(mol.GetNumAtoms(), 3)
        self.assertEqual(mol.GetNumHeavyAtoms(), 2)
        self.assertEqual(mol.GetNumBonds(), mol.nedges())
        self.assertEqual(mol.GetAtomWithIdx(1).atomic_symbol(), "O")
        self.assertEqual(mol.GetBondWithIdx(0).GetBeginAtomIdx(), 0)
        self.assertEqual([atom.atomic_symbol() for atom in mol.GetAtoms()], ["H", "O", "C"])
        self.assertEqual([(bond.a1(), bond.a2()) for bond in mol.GetBonds()], [(0, 1), (1, 2)])

        oxygen = mol.GetAtomWithIdx(1)
        mol.set_atomic_number(1, 7)
        self.assertEqual(oxygen.atomic_symbol(), "N")

    def test_molecule_common_aliases_and_counts(self):
        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        self.assertTrue(mol.ok())
        self.assertFalse(mol.empty())
        self.assertEqual(mol.GetNumAtoms(), 3)
        self.assertEqual(mol.natoms(6), 2)
        self.assertEqual(mol.natoms("O"), 1)
        with self.assertRaises(Exception):
            mol.natoms("Qq")
        self.assertEqual(mol.atomic_symbol(2), "O")
        self.assertEqual(mol.connections(1), [0, 2])
        self.assertEqual(mol.other_atom(1, 0), 0)
        self.assertEqual(mol.other_atom(1, 1), 2)

    def test_programmatic_molecule_building_and_bond_lookup(self):
        mol = lillymol_nb.Molecule()
        self.assertTrue(mol.empty())
        c1 = mol.add_atom(6)
        c2 = mol.add_atom(6)
        o = mol.add_atom(8)
        self.assertEqual((c1, c2, o), (0, 1, 2))
        self.assertEqual(mol.natoms(), 3)
        self.assertEqual(mol.add_bond(0, 1, lillymol_nb.BondType.SINGLE_BOND), 1)
        self.assertEqual(mol.add_bond(1, 2, lillymol_nb.BondType.DOUBLE_BOND), 1)
        self.assertEqual(mol.nedges(), 2)
        self.assertIsNotNone(mol.bond_between_atoms(0, 1))
        self.assertIsNone(mol.bond_between_atoms(0, 2))
        self.assertEqual(mol.bond_type_between_atoms(1, 2), lillymol_nb.BondType.DOUBLE_BOND)
        self.assertIsNone(mol.bond_type_between_atoms(0, 2))
        mol.assign_bond_numbers_to_bonds()
        self.assertTrue(mol.bond(0).bond_number_assigned())
        self.assertEqual(mol.set_bond_type_between_atoms(1, 2, lillymol_nb.BondType.SINGLE_BOND), 1)
        self.assertEqual(mol.bond_type_between_atoms(1, 2), lillymol_nb.BondType.SINGLE_BOND)
        self.assertEqual(mol.remove_bond_between_atoms(1, 2), 1)
        self.assertIsNone(mol.bond_between_atoms(1, 2))

    def test_atom_removal_and_resize(self):
        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        self.assertEqual(mol.remove_atom(2), 1)
        self.assertEqual(mol.natoms(), 2)
        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        self.assertEqual(mol.remove_atoms([0, 1, 0], 1), 1)
        self.assertEqual(mol.natoms(), 2)
        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        self.assertEqual(mol.remove_atoms(lillymol_nb.Set_of_Atoms([0, 2])), 2)
        self.assertEqual(mol.natoms(), 1)
    def test_molecule_convenience_methods(self):
        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        self.assertTrue(mol.organic_only())
        self.assertEqual(mol.non_organic_atom_count(), 0)
        self.assertTrue(mol.is_organic(0))
        self.assertTrue(mol.are_bonded(0, 1))
        self.assertFalse(mol.are_bonded(0, 2))
        self.assertEqual(mol.bonds_between(0, 2), 2)
        self.assertEqual(mol.longest_path(), 2)
        self.assertEqual(mol.most_distant_pair(), (0, 2))
        self.assertEqual(mol.atoms_on_shortest_path(0, 2), [1])
        self.assertEqual(mol.all_atoms_between(0, 2), [1])
        self.assertEqual(mol.down_the_bond(0, 1), [2])
        self.assertEqual(mol.lipinski_num_h_donors(), 1)
        self.assertEqual(mol.lipinski_num_h_acceptors(), 1)
        self.assertEqual(mol.rdkit_num_h_donors(), 1)
        self.assertEqual(mol.rdkit_num_h_acceptors(), 1)
        self.assertTrue(mol.saturated(0))
        self.assertEqual(mol.unsaturation(0), 0)

        smiles = mol.random_smiles()
        self.assertIsInstance(smiles, str)
        self.assertGreater(len(smiles), 0)
        self.assertEqual(mol.unique_kekule_smiles(), "OCC")
        self.assertEqual(mol.smiles_starting_with_atom(2), "OCC")
        order = mol.smiles_atom_order()
        self.assertEqual(sorted(order), [0, 1, 2])

        self.assertEqual(mol.renumber_atoms([2, 1, 0]), 1)
        self.assertEqual(mol.atomic_symbol(0), "O")
        with self.assertRaises(Exception):
            mol.renumber_atoms([0, 0, 1])

    def test_partial_charge_helpers(self):
        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        charges = mol.gasteiger_partial_charges()
        self.assertEqual(len(charges), mol.natoms())
        self.assertEqual(mol.partial_charge_type(), "GASTEIGER")
        self.assertAlmostEqual(mol.partial_charge(0), charges[0])
        self.assertNotEqual(sum(abs(charge) for charge in charges), 0.0)
        mol.invalidate_partial_charges()
        self.assertEqual(mol.partial_charge_type(), "")
        self.assertGreaterEqual(mol.compute_Gasteiger_partial_charges(), 0)

    def test_atoms_by_radius_single_starting_atom(self):
        mol = lillymol_nb.MolFromSmiles("CCCCC pentane")
        shells = mol.atoms_by_radius(lillymol_nb.Set_of_Atoms([2]), 3)
        self.assertEqual(len(shells), 4)
        self.assertCountEqual(shells[0], [2])
        self.assertCountEqual(shells[1], [1, 3])
        self.assertCountEqual(shells[2], [0, 4])
        self.assertEqual(shells[3], [])

    def test_atoms_by_radius_multiple_starting_atoms(self):
        mol = lillymol_nb.MolFromSmiles("CCCCC pentane")
        shells = mol.atoms_by_radius(lillymol_nb.Set_of_Atoms([0, 4]), 3)
        self.assertCountEqual(shells[0], [0, 4])
        self.assertCountEqual(shells[1], [1, 3])
        self.assertCountEqual(shells[2], [2])
        self.assertEqual(shells[3], [])

    def test_atoms_by_radius_validation(self):
        mol = lillymol_nb.MolFromSmiles("CCC propane")
        with self.assertRaises(Exception):
            mol.atoms_by_radius(lillymol_nb.Set_of_Atoms([3]), 1)
        with self.assertRaises(Exception):
            mol.atoms_by_radius(lillymol_nb.Set_of_Atoms([0]), -1)

    def test_atom_map_number_helpers(self):
        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        self.assertEqual(mol.atom_map_number(1), 0)
        mol.set_atom_map_number(1, 17)
        self.assertEqual(mol.atom_map_number(1), 17)
        self.assertEqual(mol.atom_with_atom_map_number(17), 1)
        self.assertEqual(mol.atom_with_atom_map_number(99), -1)
        mol.reset_atom_map_numbers()
        self.assertEqual(mol.atom_map_number(1), 0)

    def test_molecule_repr_debug_and_addition(self):
        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        self.assertIn("ethanol", repr(mol))
        self.assertIn("3 atoms", repr(mol))
        self.assertEqual(str(mol), "CCO ethanol")
        self.assertIn("Molecule", mol.debug_string())

        methane = lillymol_nb.MolFromSmiles("C methane")
        combined = mol + methane
        self.assertEqual(mol.natoms(), 3)
        self.assertEqual(combined.natoms(), 4)
        mol += methane
        self.assertEqual(mol.natoms(), 4)

    def test_coordinates_object_and_transforms(self):
        coords = lillymol_nb.Coordinates(3.0, 4.0, 0.0)
        self.assertAlmostEqual(coords.x(), 3.0)
        self.assertAlmostEqual(coords.y(), 4.0)
        self.assertAlmostEqual(coords.z(), 0.0)
        self.assertAlmostEqual(coords.norm(), 5.0)
        self.assertEqual(str(coords), "(3,4,0)")
        coords.normalise()
        self.assertAlmostEqual(coords.norm(), 1.0)
        coords.setxyz(1.0, 2.0, 3.0)
        coords.set_z(4.0)
        self.assertAlmostEqual(coords.length(), math.sqrt(21.0), delta=1.0e-6)
        self.assertAlmostEqual(coords.distance(lillymol_nb.Coordinates(1.0, 2.0, 4.0)), 0.0)
        self.assertAlmostEqual(coords.dot_product(lillymol_nb.Coordinates(0.0, 1.0, 0.0)), 2.0)

        mol = lillymol_nb.MolFromSmiles("CCC propane")
        mol.set_coordinates([0.0, 0.0, 0.0,
                             1.0, 0.0, 0.0,
                             2.0, 0.0, 0.0])
        mol.translate(1.0, 2.0, 3.0)
        self.assertEqual(mol.get_coordinates(), [1.0, 2.0, 3.0,
                                                 2.0, 2.0, 3.0,
                                                 3.0, 2.0, 3.0])
        mol.translate([0, 1, 0], 1, 0.0, 1.0, 0.0)
        self.assertEqual(mol.get_coordinates()[3:6], [2.0, 3.0, 3.0])
        with self.assertRaises(Exception):
            mol.translate([1, 1], 1, 0.0, 0.0, 1.0)

        ethane = lillymol_nb.MolFromSmiles("CC ethane")
        ethane.set_coordinates([0.0, 0.0, 0.0,
                                1.0, 0.0, 0.0])
        ethane.rotate(lillymol_nb.Coordinates(0.0, 0.0, 1.0), math.pi / 2.0)
        self.assertAlmostEqual(ethane.x(1), 0.0, delta=1.0e-5)
        self.assertAlmostEqual(ethane.y(1), 1.0, delta=1.0e-5)

    def test_coordinate_and_geometry_helpers(self):
        mol = lillymol_nb.MolFromSmiles("CCCO propanol")
        mol.setxyz(0, 0.0, 0.0, 0.0)
        mol.setxyz(1, 1.0, 0.0, 0.0)
        mol.setxyz(2, 1.0, 1.0, 0.0)
        mol.setxyz(3, 1.0, 1.0, 1.0)

        self.assertAlmostEqual(mol.x(1), 1.0)
        self.assertAlmostEqual(mol.y(2), 1.0)
        self.assertAlmostEqual(mol.z(3), 1.0)
        mol.setx(0, 0.25)
        mol.sety(0, 0.5)
        mol.setz(0, 0.75)
        self.assertEqual(mol.get_coordinates()[:3], [0.25, 0.5, 0.75])

        coords = [0.0, 0.0, 0.0,
                  1.0, 0.0, 0.0,
                  1.0, 1.0, 0.0,
                  1.0, 1.0, 1.0]
        mol.set_coordinates(coords)
        self.assertEqual(mol.get_coordinates(), coords)
        with self.assertRaises(Exception):
            mol.set_coordinates(coords[:-1])

        self.assertAlmostEqual(mol.distance_between_atoms(0, 1), 1.0)
        self.assertAlmostEqual(mol.bond_length(0, 1), 1.0)
        self.assertIsNone(mol.bond_length(0, 2))
        self.assertGreater(mol.bond_angle(1, 0, 2), 0.0)
        self.assertGreater(mol.dihedral_angle(0, 1, 2, 3), 0.0)
        self.assertNotEqual(mol.signed_dihedral_angle(0, 1, 2, 3), 0.0)
        self.assertAlmostEqual(mol.longest_intra_molecular_distance(), 3 ** 0.5)
        self.assertEqual(mol.highest_coordinate_dimensionality(), 3)
        self.assertGreater(mol.bump_check(2.0), 0)

    def test_canonical_and_symmetry_helpers(self):
        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        ranks = mol.canonical_ranks()
        self.assertEqual(len(ranks), mol.natoms())
        self.assertEqual(mol.canonical_rank(0), ranks[0])
        classes = [mol.symmetry_class(i) for i in range(mol.natoms())]
        self.assertEqual(mol.number_symmetry_classes(), len(set(classes)))

        ethane = lillymol_nb.MolFromSmiles("CC ethane")
        equivalents = ethane.symmetry_equivalents(0)
        self.assertIn(1, equivalents)

    def test_remove_helpers(self):
        mol = lillymol_nb.MolFromSmiles("[Na]OC sodium_ethoxide")
        self.assertEqual(mol.non_organic_atom_count(), 1)
        self.assertEqual(mol.remove_non_periodic_table_elements(), 0)

        mol = lillymol_nb.MolFromSmiles("CC ethane")
        self.assertEqual(mol.AddHs(), 6)
        self.assertEqual(mol.remove_explicit_hydrogens(), 2)
        self.assertEqual(mol.natoms(), 2)

        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        self.assertTrue(mol.remove_bonds_to_atom(1))
        self.assertEqual(mol.nedges(), 0)

        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        self.assertEqual(mol.remove_edge(0), 1)
        self.assertEqual(mol.nedges(), 1)
        self.assertEqual(mol.chop(1), 2)
        self.assertEqual(mol.natoms(), 2)

    def test_charges_isotopes_and_hydrogens(self):
        mol = lillymol_nb.MolFromSmiles("C[NH3+] methylammonium")
        self.assertTrue(mol.has_formal_charges())
        self.assertEqual(mol.number_formal_charges(), 1)
        self.assertEqual(mol.net_formal_charge(), 1)
        self.assertEqual(mol.formal_charge(1), 1)
        mol.set_formal_charge(1, 0)
        self.assertEqual(mol.formal_charge(1), 0)
        self.assertEqual(mol.net_formal_charge(), 0)

        mol = lillymol_nb.MolFromSmiles("[3CH3]CO labelled")
        self.assertEqual(mol.number_isotopic_atoms(), 1)
        self.assertEqual(mol.first_atom_with_isotope(3), 0)
        self.assertEqual(mol.transform_to_non_isotopic_form(), 1)
        self.assertEqual(mol.number_isotopic_atoms(), 0)

        mol = lillymol_nb.MolFromSmiles("CC ethane")
        self.assertEqual(mol.implicit_hydrogens(0), 3)
        self.assertEqual(mol.explicit_hydrogens(0), 0)
        self.assertEqual(mol.hcount(0), 3)
        self.assertTrue(mol.valence_ok())
        self.assertTrue(mol.valence_ok(0))
        self.assertEqual(mol.AddHs(), 6)
        self.assertEqual(mol.natoms(), 8)
        self.assertGreater(mol.RemoveHs(), 0)
        self.assertEqual(mol.natoms(), 2)

    def test_ring_aliases_and_aromaticity(self):
        mol = lillymol_nb.MolFromSmiles("c1ccccc1 benzene")
        self.assertTrue(mol.IsInRing(0))
        self.assertTrue(mol.in_ring_of_size(0, 6))
        self.assertTrue(mol.IsAtomInRingOfSize(0, 6))
        self.assertEqual(mol.NumAtomRings(0), 1)
        self.assertTrue(mol.is_aromatic(0))
        self.assertEqual(mol.aromatic_atom_count(), 6)
        self.assertEqual(mol.aromatic_ring_count(), 1)

    def test_chiral_centre_access_and_helpers(self):
        mol = lillymol_nb.MolFromSmiles("F[C@H](Cl)Br chiral")
        self.assertEqual(mol.number_chiral_centres(), 1)
        centre = mol.chiral_centre_at_atom(1)
        self.assertIsNotNone(centre)
        self.assertEqual(centre.atom(), 1)
        self.assertTrue(centre.involves(0))
        self.assertEqual(centre.implicit_hydrogens(), 1)
        self.assertEqual(centre.lone_pairs(), 0)
        self.assertIn("<Chiral_Centre atom 1", repr(centre))
        centres = mol.chiral_centres()
        self.assertEqual(len(centres), 1)
        self.assertEqual(centres[0].atom(), 1)
        self.assertTrue(lillymol_nb.is_actually_chiral(mol, 1))
        self.assertFalse(lillymol_nb.is_actually_chiral(mol, 0))
        tag = lillymol_nb.tetrahedral_chirality(mol, 1)
        self.assertIn(tag, [lillymol_nb.ChiralType.CHI_TETRAHEDRAL_CW,
                            lillymol_nb.ChiralType.CHI_TETRAHEDRAL_CCW])
        self.assertEqual(lillymol_nb.tetrahedral_chirality(mol, 0), None)
        self.assertEqual(lillymol_nb.tetrahedral_chirality(mol, 1, check_is_chiral=True), tag)
        self.assertEqual(mol.invert_chirality_on_atom(1), 1)
        inverted = lillymol_nb.tetrahedral_chirality(mol, 1)
        self.assertIn(inverted, [lillymol_nb.ChiralType.CHI_TETRAHEDRAL_CW,
                                 lillymol_nb.ChiralType.CHI_TETRAHEDRAL_CCW])
        self.assertNotEqual(inverted, tag)
        self.assertEqual(mol.remove_chiral_centre_at_atom(1), 1)
        self.assertEqual(mol.number_chiral_centres(), 0)
        self.assertIsNone(mol.chiral_centre_at_atom(1))

    def test_remove_all_chiral_centres(self):
        mol = lillymol_nb.MolFromSmiles("F[C@H](Cl)Br chiral")
        self.assertEqual(mol.number_chiral_centres(), 1)
        self.assertEqual(mol.remove_all_chiral_centres(), 1)
        self.assertEqual(mol.number_chiral_centres(), 0)
        self.assertIsNone(lillymol_nb.tetrahedral_chirality(mol, 1))

    def test_atom_scalars_and_vectors(self):
        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        self.assertEqual(mol.atomic_number(0), 6)
        self.assertEqual(mol.atomic_numbers(), [6, 6, 8])
        self.assertEqual(mol.ncon(1), 2)
        self.assertEqual(mol.nbonds(1), 2)
        self.assertEqual(mol.isotopes(), [0, 0, 0])
        self.assertEqual(mol.set_isotope(1, 7), 1)
        self.assertEqual(mol.isotope(1), 7)
        self.assertEqual(mol.isotopes(), [0, 7, 0])
        self.assertEqual(mol.remove_isotopes(), 1)
        self.assertEqual(mol.isotopes(), [0, 0, 0])

    def test_fragments_and_formula(self):
        mol = lillymol_nb.MolFromSmiles("CC.O mixture")
        self.assertEqual(mol.number_fragments(), 2)
        self.assertEqual(mol.atoms_in_largest_fragment(), 2)
        self.assertEqual(mol.molecular_formula(), "C2OH8")

    def test_molecular_weight(self):
        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        self.assertAlmostEqual(mol.amw(), 46.069, delta=0.001)
        self.assertAlmostEqual(lillymol_nb.molecular_weight(mol), 46.069, delta=0.001)
        self.assertGreater(mol.exact_mass(), 0.0)

    def test_atom_access(self):
        atom = lillymol_nb.Atom(6)
        self.assertEqual(atom.atomic_number(), 6)
        self.assertEqual(atom.atomic_symbol(), "C")
        self.assertEqual(atom.ncon(), 0)
        self.assertEqual(atom.nbonds(), 0)
        self.assertEqual(atom.formal_charge(), 0)
        self.assertTrue(atom.is_organic())
        self.assertGreater(atom.atomic_weight(), 12.0)
        self.assertEqual(atom.exact_mass(), 12.0)
        self.assertIn("<Atom C", repr(atom))

        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        carbon = mol.atom(1)
        self.assertEqual(carbon.atomic_number(), 6)
        self.assertEqual(carbon.atomic_symbol(), "C")
        self.assertEqual(len(carbon), 2)
        self.assertEqual(carbon.connections(1), [0, 2])
        self.assertEqual([(bond.a1(), bond.a2()) for bond in carbon], [(0, 1), (1, 2)])
        self.assertEqual(carbon[0].other(1), 0)
        self.assertTrue(carbon.is_bonded_to(0))
        self.assertIn(2, carbon)
        self.assertEqual(carbon.other(1, 0), 0)
        self.assertTrue(carbon.valence_ok())
        self.assertEqual(carbon.unsaturation(), 0)
        self.assertTrue(carbon.fully_saturated())

    def test_atom_coordinate_distance_helpers(self):
        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        mol.set_coordinates([0.0, 0.0, 0.0,
                             1.0, 0.0, 0.0,
                             1.0, 1.0, 0.0])
        a0 = mol.atom(0)
        a1 = mol.atom(1)
        self.assertAlmostEqual(a0.x(), 0.0)
        self.assertAlmostEqual(a1.x(), 1.0)
        self.assertAlmostEqual(a0.distance(a1), 1.0)
        self.assertAlmostEqual(a0.distance(lillymol_nb.Coordinates(0.0, 1.0, 0.0)), 1.0)
        self.assertAlmostEqual(a0 - a1, 1.0)
        self.assertIn("<Atom C", str(a0))

    def test_atom_view_reflects_parent_mutation(self):
        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        atom = mol.atom(1)
        self.assertEqual(atom.atomic_number(), 6)
        self.assertEqual(atom.atomic_symbol(), "C")
        self.assertEqual(mol.set_atomic_number(1, 7), 1)
        self.assertEqual(atom.atomic_number(), 7)
        self.assertEqual(atom.atomic_symbol(), "N")
        self.assertEqual(mol.atomic_number(1), 7)

    def test_bond_access(self):
        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        self.assertEqual(mol.nedges(), 2)
        bond = mol.bond(0)
        self.assertEqual(bond.a1(), 0)
        self.assertEqual(bond.a2(), 1)
        self.assertEqual(bond.other(0), 1)
        self.assertTrue(bond.involves(1))
        self.assertTrue(bond.is_single_bond())
        self.assertFalse(bond.is_double_bond())
        self.assertFalse(bond.is_triple_bond())
        self.assertFalse(bond.is_aromatic())
        self.assertIn("<Bond 0-1>", repr(bond))

        bonds = mol.bonds()
        self.assertEqual(len(bonds), 2)
        self.assertEqual([(b.a1(), b.a2()) for b in bonds], [(0, 1), (1, 2)])

    def test_bond_view_reflects_parent_mutation(self):
        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        bond = mol.bond(0)
        self.assertTrue(bond.is_single_bond())
        self.assertEqual(bond.btype(), lillymol_nb.BondType.SINGLE_BOND)
        self.assertEqual(mol.set_bond_type_between_atoms(0, 1, lillymol_nb.BondType.DOUBLE_BOND), 1)
        self.assertFalse(bond.is_single_bond())
        self.assertTrue(bond.is_double_bond())
        self.assertEqual(bond.btype(), lillymol_nb.BondType.DOUBLE_BOND)

    def test_bond_rdkit_style_aliases_and_contains(self):
        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        bond = mol.bond(0)
        self.assertEqual(bond.GetBeginAtomIdx(), 0)
        self.assertEqual(bond.GetEndAtomIdx(), 1)
        self.assertEqual(bond.GetBondType(), lillymol_nb.BondType.SINGLE_BOND)
        self.assertIn(0, bond)
        self.assertNotIn(2, bond)
        self.assertIn("<Bond 0-1>", str(bond))

    def test_bond_ring_membership(self):
        mol = lillymol_nb.MolFromSmiles("c1ccccc1 benzene")
        bond = mol.bond(0)
        mol.ring_membership()
        self.assertEqual(bond.nrings(), 1)
        self.assertTrue(bond.IsInRing())


    def test_set_of_atoms(self):
        atoms = lillymol_nb.Set_of_Atoms([1, 3])
        self.assertFalse(atoms.empty())
        self.assertEqual(len(atoms), 2)
        self.assertEqual(atoms.size(), 2)
        self.assertIn(1, atoms)
        self.assertNotIn(2, atoms)
        self.assertEqual(atoms[0], 1)
        self.assertEqual(list(atoms), [1, 3])
        atoms[0] = 0
        self.assertEqual(atoms.as_list(), [0, 3])
        atoms.append(4)
        atoms.extend([5, 6])
        self.assertEqual(atoms.as_list(), [0, 3, 4, 5, 6])
        self.assertTrue(atoms.contains_both(3, 5))
        self.assertEqual(atoms.atoms_in_common(lillymol_nb.Set_of_Atoms([5, 9])), 1)
        self.assertEqual(atoms.first_atom_in_common(lillymol_nb.Set_of_Atoms([5, 9])), 5)
        atoms += lillymol_nb.Set_of_Atoms([7])
        atoms += 8
        self.assertEqual(atoms.as_list()[-2:], [7, 8])
        self.assertEqual(lillymol_nb.Set_of_Atoms([1, 2]), [1, 2])

    def test_ring_info_helpers(self):
        mol = lillymol_nb.MolFromSmiles("c1ccccc1 benzene")
        ring_info = mol.GetRingInfo()
        self.assertEqual(ring_info.NumRings(), 1)
        self.assertEqual(ring_info.num_rings(), 1)
        self.assertEqual([sorted(ring) for ring in ring_info.AtomRings()], [[0, 1, 2, 3, 4, 5]])
        bond_rings = ring_info.BondRings()
        self.assertEqual(len(bond_rings), 1)
        self.assertEqual(len(bond_rings[0]), 6)
        self.assertEqual(ring_info.NumAtomRings(0), 1)
        self.assertEqual(ring_info.NumBondRings(0), 1)
        self.assertTrue(ring_info.IsAtomInRingOfSize(0, 6))
        self.assertTrue(ring_info.AreAtomsInSameRing(0, 3))
        self.assertTrue(ring_info.AreBondsInSameRing(0, 2))
        with self.assertRaises(Exception):
            ring_info.NumBondRings(99)

        fused = lillymol_nb.MolFromSmiles("c1ccc2ccccc2c1 naphthalene")
        fused_info = fused.ring_info()
        self.assertEqual(fused_info.NumRings(), 2)
        self.assertEqual(len(fused_info.AtomRings()), 2)
        self.assertEqual(len(fused_info.BondRings()), 2)
        self.assertGreaterEqual(fused_info.NumAtomRings(3), 1)

    def test_ring_access(self):
        mol = lillymol_nb.MolFromSmiles("c1ccccc1 benzene")
        self.assertEqual(mol.nrings(), 1)
        self.assertEqual(mol.get_ring_membership(), [1, 1, 1, 1, 1, 1])
        self.assertTrue(mol.is_ring_atom(0))
        self.assertEqual(mol.ring_bond_count(0), 2)
        self.assertEqual(mol.ring_bond_counts(), [2, 2, 2, 2, 2, 2])
        self.assertEqual(mol.number_ring_systems(), 1)
        self.assertEqual(mol.largest_ring_size(), 6)
        self.assertTrue(mol.in_same_ring(0, 3))
        self.assertTrue(mol.in_same_ring_system(0, 3))
        self.assertEqual(mol.fused_system_size(0), 1)

        ring = mol.ring(0)
        self.assertEqual(ring.size(), 6)
        self.assertEqual(len(ring), 6)
        self.assertEqual(ring.ring_number(), 0)
        self.assertTrue(ring.contains_bond(0, 1))
        self.assertTrue(ring.contains_both(0, 3))
        self.assertFalse(ring.is_fused())
        self.assertEqual(ring.fused_ring_neighbours(), 0)
        self.assertEqual(set(ring.as_list()), set(range(6)))
        self.assertEqual(set(ring), set(range(6)))

        rings = mol.rings()
        self.assertEqual(len(rings), 1)
        self.assertEqual(rings[0].as_list(), ring.as_list())


    def test_molecule_sequence_and_copy(self):
        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        self.assertEqual(len(mol), 3)
        self.assertEqual(mol[2].atomic_symbol(), "O")
        self.assertEqual([atom.atomic_symbol() for atom in mol], ["C", "C", "O"])
        self.assertEqual([atom.ncon() for atom in mol], [1, 2, 1])

        mol_copy = copy.copy(mol)
        self.assertEqual(mol_copy.smiles(), mol.smiles())
        mol_copy.set_name("copy")
        self.assertEqual(mol.name(), "ethanol")
        self.assertEqual(mol_copy.name(), "copy")


    def test_molecule_contains(self):
        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        self.assertIn(6, mol)
        self.assertNotIn(7, mol)
        self.assertIn("C", mol)
        self.assertIn("O", mol)
        self.assertNotIn("N", mol)
        with self.assertRaises(Exception):
            "Qq" in mol

        carbon = lillymol_nb.QueryFromSmarts("C")
        nitrogen = lillymol_nb.QueryFromSmarts("N")
        self.assertIn(carbon, mol)
        self.assertNotIn(nitrogen, mol)


    def test_substructure_query_object(self):
        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        qry = lillymol_nb.SubstructureQuery()
        self.assertTrue(qry.build_from_smarts("C"))
        self.assertEqual(qry.number_elements(), 1)
        self.assertEqual(qry.active(), 1)
        self.assertEqual(qry.substructure_search(mol), 2)
        matches = qry.substructure_search_matches(mol)
        self.assertEqual([m.as_list() for m in matches], [[0], [1]])
        self.assertEqual(qry.substructure_search_match_lists(mol), [[0], [1]])

        qry.set_max_matches_to_find(1)
        self.assertEqual(qry.substructure_search(mol), 1)

    def test_substructure_results_object(self):
        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        qry = lillymol_nb.QueryFromSmarts("C")
        self.assertIsNotNone(qry)
        sresults = lillymol_nb.SubstructureResults()
        self.assertEqual(qry.substructure_search(mol, sresults), 2)
        self.assertEqual(sresults.number_embeddings(), 2)
        self.assertEqual(len(sresults), 2)
        self.assertEqual(sresults.embeddings(), [[0], [1]])
        self.assertEqual([embedding.as_list() for embedding in sresults], [[0], [1]])
        self.assertEqual(sresults[0].as_list(), [0])
        self.assertEqual(sresults.each_embedding_set_vector(mol.natoms(), 1), [1, 1, 0])
        self.assertGreaterEqual(sresults.max_query_atoms_matched_in_search(), 1)
        sresults.reset()
        self.assertEqual(sresults.number_embeddings(), 0)

    def test_substructure_free_functions(self):
        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        self.assertTrue(lillymol_nb.HasSubstructMatch(mol, "CO"))
        self.assertFalse(lillymol_nb.HasSubstructMatch(mol, "N"))
        self.assertEqual(lillymol_nb.CountSubstructMatches(mol, "C"), 2)
        self.assertEqual(lillymol_nb.CountSubstructMatches(mol, "C", max_matches_to_find=1), 1)
        self.assertEqual(lillymol_nb.GetSubstructMatches(mol, "CO"), [[1, 2]])
        self.assertEqual(lillymol_nb.GetSubstructMatches(mol, "N"), [])
        self.assertIsNone(lillymol_nb.QueryFromSmarts("["))
        with self.assertRaises(Exception):
            lillymol_nb.HasSubstructMatch(mol, "[")


    def test_element_transformations(self):
        etrans = lillymol_nb.ElementTransformations()
        self.assertFalse(etrans.active())
        self.assertTrue(etrans.add("I=Cl"))
        self.assertTrue(etrans.add("Br=Cl"))
        self.assertTrue(etrans.active())
        self.assertFalse(etrans.add("not a transform"))

        mol = lillymol_nb.MolFromSmiles("ICBr halides")
        self.assertEqual(mol.natoms("I"), 1)
        self.assertEqual(mol.natoms("Br"), 1)
        self.assertEqual(etrans.process(mol), 2)
        self.assertEqual(mol.natoms("I"), 0)
        self.assertEqual(mol.natoms("Br"), 0)
        self.assertEqual(mol.natoms("Cl"), 2)

    def test_standardise(self):
        mol = lillymol_nb.MolFromSmiles("CC(=O)[O-] acetate")
        standardise = lillymol_nb.Standardise()
        self.assertEqual(standardise.process(mol), 0)
        self.assertEqual(mol.smiles(), "CC(=O)[O-]")

        standardise.activate_all()
        self.assertEqual(standardise.process(mol), 1)
        self.assertEqual(mol.smiles(), "CC(=O)O")

    def test_fingerprint_default_and_tanimoto(self):
        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        bits = lillymol_nb.linear_fingerprint(mol)
        self.assertIsInstance(bits, list)
        self.assertEqual(len(bits), 2048)
        self.assertGreater(sum(bits), 0)
        self.assertAlmostEqual(lillymol_nb.tanimoto(bits, bits), 1.0)

        bits2 = lillymol_nb.linear_fingerprint(mol, nbits=512, atype_specification="")
        self.assertIsNotNone(bits2)
        self.assertEqual(len(bits2), 512)
        self.assertIsNone(lillymol_nb.linear_fingerprint(mol, nbits=512, atype_specification="BAD"))

    def test_fingerprint_creators(self):
        mol = lillymol_nb.MolFromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C caffeine")
        ecfp = lillymol_nb.ECFingerprintCreator(512)
        bits = ecfp.fingerprint(mol)
        self.assertIsInstance(bits, list)
        self.assertEqual(len(bits), 512)
        self.assertGreater(sum(bits), 0)

        same_molecule = lillymol_nb.MolFromSmiles("Cn1c(=O)c2c(ncn2C)n(C)c1=O caffeine")
        self.assertEqual(mol.unique_smiles(), same_molecule.unique_smiles())
        self.assertEqual(bits, ecfp.fingerprint(same_molecule))

        linear = lillymol_nb.LinearFingerprintCreator(256)
        linear.set_max_length(5)
        self.assertEqual(len(linear.fingerprint(mol)), 256)

        atom_pair = lillymol_nb.AtomPairFingerprintCreator(256)
        atom_pair.set_min_separation(1)
        atom_pair.set_max_separation(5)
        self.assertEqual(len(atom_pair.fingerprint(mol)), 256)

    def test_descriptor_helpers(self):
        mol = lillymol_nb.MolFromSmiles("CCO ethanol")
        self.assertIsNotNone(lillymol_nb.alogp(mol))
        self.assertIsNotNone(lillymol_nb.xlogp(mol))
        self.assertIsNotNone(lillymol_nb.tpsa(mol))
        self.assertEqual(lillymol_nb.HbaHbd(mol), (1, 1))
        self.assertEqual(lillymol_nb.NumHAcceptors(mol), mol.lipinski_num_h_acceptors())
        self.assertEqual(lillymol_nb.NumHDonors(mol), mol.lipinski_num_h_donors())
        self.assertEqual(lillymol_nb.RDKitNumHAcceptors(mol), mol.rdkit_num_h_acceptors())
        self.assertEqual(lillymol_nb.RDKitNumHDonors(mol), mol.rdkit_num_h_donors())
        self.assertAlmostEqual(lillymol_nb.fraction_csp3(mol), 1.0)

        benzene = lillymol_nb.MolFromSmiles("c1ccccc1 benzene")
        self.assertAlmostEqual(lillymol_nb.fraction_csp3(benzene), 0.0)
        water = lillymol_nb.MolFromSmiles("O water")
        self.assertAlmostEqual(lillymol_nb.fraction_csp3(water), 0.0)


class TestNanobindTSubstructure(LillyMolNanobindTestCase):

    def test_no_queries(self):
        ts = lillymol_nb.TSubstructure()
        self.assertFalse(ts.substructure_search("C methane"))

    def test_single_query_single_molecule(self):
        ts = lillymol_nb.TSubstructure()
        self.assertTrue(ts.add_query_from_smarts("C carbon"))
        mol = lillymol_nb.MolFromSmiles("CC")
        self.assertTrue(ts.substructure_search(mol))
        self.assertEqual(ts.num_matches(mol), [2])
        self.assertEqual(ts.number_queries(), 1)
        self.assertEqual(ts.query_names(), ["carbon"])

    def test_multiple_query_single_molecule(self):
        ts = lillymol_nb.TSubstructure()
        self.assertTrue(ts.add_query_from_smarts("C carbon"))
        self.assertTrue(ts.add_query_from_smarts("N nitrogen"))
        mol = lillymol_nb.MolFromSmiles("CC")
        self.assertEqual(ts.num_matches(mol), [2, 0])

    def test_batch_substructure_search_molecules(self):
        ts = lillymol_nb.TSubstructure()
        self.assertTrue(ts.add_query_from_smarts("C carbon"))
        smiles = ["C methane", "CC ethane", "CCC propane", "C1CC1 cyclopropane", "c1ccccc1 benzene"]
        mols = [lillymol_nb.MolFromSmiles(smi) for smi in smiles]
        self.assertEqual(ts.substructure_search(mols), [True, True, True, True, False])

    def test_batch_substructure_search_smiles(self):
        ts = lillymol_nb.TSubstructure()
        self.assertTrue(ts.add_query_from_smarts("N nitrogen"))
        smiles = ["C methane", "N nitrogen", "CN methylamine"]
        self.assertEqual(ts.substructure_search(smiles), [False, True, True])

    def test_batch_num_matches_molecules(self):
        ts = lillymol_nb.TSubstructure()
        self.assertTrue(ts.add_query_from_smarts("C carbon"))
        self.assertTrue(ts.add_query_from_smarts("N nitrogen"))
        smiles = ["C methane", "CC ethane", "N nitrogen", "O oxygen", "CN CN"]
        mols = [lillymol_nb.MolFromSmiles(smi) for smi in smiles]
        self.assertEqual(ts.num_matches(mols), [[1, 0], [2, 0], [0, 1], [0, 0], [1, 1]])

    def test_must_match_all_queries(self):
        ts = lillymol_nb.TSubstructure()
        self.assertTrue(ts.add_query_from_smarts("C carbon"))
        self.assertTrue(ts.add_query_from_smarts("N nitrogen"))
        mol = lillymol_nb.MolFromSmiles("C")
        self.assertTrue(ts.substructure_search(mol))
        ts.must_match_all_queries = True
        self.assertFalse(ts.substructure_search(mol))

    def test_unique_embeddings_only(self):
        ts = lillymol_nb.TSubstructure()
        mol = lillymol_nb.MolFromSmiles("CC(C)(C)c1c(C(C)(C)C)c(C(C)(C)C)c(C(C)(C)C)c(C(C)(C)C)c1C(C)(C)C")
        self.assertTrue(ts.add_query_from_smarts("CC(C)(C)c1c(C(C)(C)C)c(C(C)(C)C)c(C(C)(C)C)c(C(C)(C)C)c1C(C)(C)C"))
        self.assertEqual(ts.num_matches(mol), [559872])
        ts.set_unique_embeddings_only(True)
        self.assertEqual(ts.num_matches(mol), [1])

    def test_max_matches_to_find(self):
        ts = lillymol_nb.TSubstructure()
        mol = lillymol_nb.MolFromSmiles("c1ccccc1")
        self.assertTrue(ts.add_query_from_smarts("c1ccccc1"))
        self.assertEqual(ts.num_matches(mol), [12])
        ts.set_max_matches_to_find(5)
        self.assertEqual(ts.num_matches(mol), [5])

    def test_reduce_to_largest_fragment(self):
        ts = lillymol_nb.TSubstructure()
        mol = lillymol_nb.MolFromSmiles("CC.C")
        self.assertTrue(ts.add_query_from_smarts("C"))
        self.assertEqual(ts.num_matches(mol), [3])
        ts.set_reduce_to_largest_fragment(True)
        self.assertEqual(ts.num_matches(mol), [2])

    def test_make_implicit_hydrogens_explicit(self):
        ts = lillymol_nb.TSubstructure()
        mol = lillymol_nb.MolFromSmiles("CC")
        self.assertTrue(ts.add_query_from_smarts("[#1]-C"))
        self.assertEqual(ts.num_matches(mol), [0])
        ts.set_make_implicit_hydrogens_explicit(True)
        self.assertEqual(ts.num_matches(mol), [6])

    def test_label_matched_atoms(self):
        ts = lillymol_nb.TSubstructure()
        mol = lillymol_nb.MolFromSmiles("Cc1ccccc1")
        self.assertTrue(ts.add_query_from_smarts("C"))
        self.assertTrue(ts.add_query_from_smarts("c"))
        ts.isotope = 4
        self.assertEqual(ts.label_matched_atoms(mol), 2)
        self.assertEqual(mol.unique_smiles(), "[4CH3][4c]1[4cH][4cH][4cH][4cH][4cH]1")

    def test_matched_atoms(self):
        ts = lillymol_nb.TSubstructure()
        self.assertTrue(ts.add_query_from_smarts("C"))
        self.assertTrue(ts.add_query_from_smarts("N"))
        mol = lillymol_nb.MolFromSmiles("CC")
        self.assertEqual(ts.matched_atoms(mol), [[0, 1], []])


if __name__ == "__main__":
    unittest.main()
