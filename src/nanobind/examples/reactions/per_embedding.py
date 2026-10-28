"""Generate separate products from the two chlorine sites of one scaffold."""
from lillymol import MolFromSmiles, Reaction


def main():
    reaction = Reaction()
    assert reaction.construct_from_smirks("[Cl:1]>>[F:1]")
    reaction.set_name("One chlorine to fluorine")
    scaffold = MolFromSmiles("ClCCC(Cl)C dichloropentane")
    original = scaffold.unique_smiles()
    matches = reaction.substructure_search_matches(scaffold)
    assert matches is not None
    assert len(matches) == 2
    products = []
    for match in matches:
        product = reaction.perform_reaction(scaffold, match)
        assert product is not None
        product.set_name(scaffold.name() + " " + reaction.name())
        products.append(product.unique_smiles())
        print(product.smiles(), product.name())
    expected = {MolFromSmiles(s).unique_smiles()
                for s in ("FCCC(Cl)C", "ClCCC(F)C")}
    assert set(products) == expected
    assert scaffold.unique_smiles() == original
    assert reaction.substructure_search_matches(MolFromSmiles("CCO")) is None


if __name__ == "__main__":
    main()
