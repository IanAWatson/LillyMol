"""Enumerate two acids against two stored amines; skip a nonmatching scaffold."""
from pathlib import Path

from lillymol import (MolFromSmiles, Reaction, ReactionIterator,
                     SidechainMatchConditions)


def main():
    directory = Path(__file__).parent
    reaction = Reaction()
    assert reaction.read(str(directory / "amide.textproto"))
    conditions = SidechainMatchConditions()
    count = reaction.add_sidechain_reagents(
        0, str(directory / "amines.smi"), conditions)
    assert count == 2
    assert reaction.sidechain_name(0, 0) == "methylamine"
    assert reaction.number_sidechains_with_reagents() == 1
    iterator = ReactionIterator(reaction)
    products = []
    for smiles in ("CC(=O)O acetic_acid", "CCC(=O)O propionic_acid", "CCC no_acid"):
        scaffold = MolFromSmiles(smiles)
        matches = reaction.substructure_search_matches(scaffold)
        if matches is None:
            continue
        iterator.reset()
        while iterator.active():
            names = reaction.reagent_names(iterator)
            for match in matches:
                product = reaction.perform_reaction(scaffold, match, iterator)
                assert product is not None
                products.append(product.unique_smiles())
                print(product.smiles(), scaffold.name(), *names)
            iterator.increment()
    expected = {MolFromSmiles(s).unique_smiles() for s in
                ("CC(=O)NC", "CC(=O)NCC", "CCC(=O)NC", "CCC(=O)NCC")}
    assert len(products) == 4
    assert set(products) == expected
    assert reaction.remove_all_reagents() == 1  # Number of sidechains, not reagents.
    assert reaction.number_sidechains_with_reagents() == 0


if __name__ == "__main__":
    main()
