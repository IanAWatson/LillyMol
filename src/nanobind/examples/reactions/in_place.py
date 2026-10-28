"""Replace both chlorine atoms in one molecule; check the resulting structure."""
from pathlib import Path

from lillymol import MolFromSmiles, Reaction


def main():
    reaction = Reaction()
    reaction_file = Path(__file__).with_name("chlorine_to_fluorine.textproto")
    assert reaction.read(str(reaction_file))
    molecule = MolFromSmiles("ClCCCCl dichloropropane")
    assert reaction.in_place_transformations(molecule)
    assert molecule.unique_smiles() == MolFromSmiles("FCCCF").unique_smiles()
    assert not reaction.in_place_transformations(molecule)  # No chlorine remains.
    print(reaction)
    print(molecule.smiles(), molecule.name())


if __name__ == "__main__":
    main()
