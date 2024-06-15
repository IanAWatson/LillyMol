from lillymol import *

mol = MolFromSmiles("CCOC(=O)C(C)(C)OC1=CC=C(C=C1)Cl.CO.C1=CC(=CC=C1C(=O)N[C@@H](CCC(=O)O)C(=O)O)NCC2=CN=C3C(=N2)C(=O)NC(=N3)N")
mol.reduce_to_largest_fragment_carefully()
print(mol.aromatic_smiles())

mol = MolFromSmiles("CCOC(=O)C(C)(C)OC1=CC=C(C=C1)Cl.CO.C1=CC(=CC=C1C(=O)N[C@@H](CCC(=O)O)C(=O)O)NCC2=CN=C3C(=N2)C(=O)NC(=N3)N")
components = mol.create_components()
largest = max(components, key=lambda m: m.natoms())
print(largest.aromatic_smiles())
