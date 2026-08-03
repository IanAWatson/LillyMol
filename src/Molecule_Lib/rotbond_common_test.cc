// rotbond_common.h does not include molecule.h, it relies on the caller having
// done so already.
#include "Molecule_Lib/molecule.h"

#include "Molecule_Lib/rotbond_common.h"

#include "gtest/gtest.h"

namespace {

struct SmilesRotbond {
  const char* smiles;
  int rotbond;
};

class TestExpensiveRotbond : public testing::TestWithParam<SmilesRotbond> {
  protected:
    Molecule _m;
    quick_rotbond::QuickRotatableBonds _rotbond;

    void SetUp() override {
      _rotbond.set_calculation_type(quick_rotbond::QuickRotatableBonds::RotBond::kExpensive);
    }
};

TEST_P(TestExpensiveRotbond, Count) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));

  EXPECT_EQ(_rotbond.Process(_m), params.rotbond) << params.smiles;
}

// The conjugated linkage exclusions. A single bond from a carbonyl like centre
// to a heteroatom is not freely rotatable. These used to require the partner to
// be Nitrogen, so the ester, thioester and amidine cases counted as rotatable
// and the totals ran above RDKit's.
INSTANTIATE_TEST_SUITE_P(ConjugatedLinkages, TestExpensiveRotbond, testing::Values(
  // A plain chain and a plain ether, for contrast - nothing excluded.
  SmilesRotbond{"CCCC", 1},
  SmilesRotbond{"CCOCC", 2},

  // Amide and thioamide. These were always excluded.
  SmilesRotbond{"CC(=O)NCC", 1},
  SmilesRotbond{"CC(=S)NCC", 1},

  // Ester, thioester, amidine and carbamate. These are the ones that changed.
  SmilesRotbond{"CC(=O)OCC", 1},
  SmilesRotbond{"CC(=O)SCC", 1},
  SmilesRotbond{"CC(=N)NCC", 1},
  SmilesRotbond{"CCOC(=O)NCC", 2},

  // A sulfonamide S-N bond is excluded. RDKit counts it, because its pattern
  // requires the doubly bonded centre to be carbon. Deliberate difference.
  SmilesRotbond{"CS(=O)(=O)NCC", 1},

  // Aromatic rings are held in a Kekule form, so an aromatic carbon really does
  // have a double bond to an aromatic ring Nitrogen. That must not make the
  // exocyclic bond look conjugated, or every aminopyridine would be wrong.
  SmilesRotbond{"CCNc1ccccn1", 2},
  SmilesRotbond{"CCNc1ccccc1", 2},

  // A lactam C=O is exocyclic even though the carbon is in the ring, so the
  // linkage is still recognised.
  SmilesRotbond{"CCNC(=O)C1CCCC1", 2},

  // The other exclusions, unchanged - terminal atoms, triple bonds, ring bonds,
  // trihalomethyl and t-butyl.
  SmilesRotbond{"CCC#N", 0},
  SmilesRotbond{"C1CCCCC1", 0},
  SmilesRotbond{"FC(F)(F)CC", 0},
  SmilesRotbond{"CC(C)(C)CC", 0}
));

}  // namespace
