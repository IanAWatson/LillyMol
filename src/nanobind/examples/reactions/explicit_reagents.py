"""Form one amide from an explicitly supplied scaffold and sidechain reagent."""
from pathlib import Path

from lillymol import MolFromSmiles, Reaction


def main():
    reaction = Reaction()
    text = Path(__file__).with_name("amide.textproto").read_text()
    assert reaction.construct_from_textproto(text)
    assert reaction.number_sidechains() == 1
    scaffold = MolFromSmiles("CC(=O)O acetic_acid")
    amine = MolFromSmiles("CN methylamine")
    original = [m.unique_smiles() for m in (scaffold, amine)]
    product = reaction.perform_reaction([scaffold, amine])
    assert product is not None
    assert product.unique_smiles() == MolFromSmiles("CC(=O)NC").unique_smiles()
    assert [m.unique_smiles() for m in (scaffold, amine)] == original
    assert reaction.number_sidechains_with_reagents() == 0
    print(product.smiles())


if __name__ == "__main__":
    main()
