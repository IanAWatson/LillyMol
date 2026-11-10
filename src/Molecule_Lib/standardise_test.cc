#include <algorithm>
#include <utility>
#include <vector>

#include "googletest/include/gtest/gtest.h"

#include "aromatic.h"
#include "molecule.h"
#include "smiles.h"
#include "standardise.h"

namespace {


class TestStandardisation : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    Chemical_Standardisation _chemical_standardisation;

    IWString _smiles;
    Molecule _m1;
    Molecule _m2;
};

void
TestStandardisation::SetUp()
{
  set_global_aromaticity_type(Daylight);
  set_unique_smiles_legacy_atom_ordering(true);
}

TEST_F(TestStandardisation, EmptyMolecule)
{
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.activate_all();
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestAcidYes)
{
  _smiles = "CC(=O)[O-]";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  EXPECT_EQ(_m1.molecular_formula(), "C2O2H3");
  EXPECT_EQ(_m1.smiles(), "CC(=O)[O-]");

  _chemical_standardisation.Activate(CS_ACID, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.molecular_formula(), "C2O2H4");
  EXPECT_EQ(_m1.smiles(), "CC(=O)O");
  EXPECT_EQ(_m1.unique_smiles(), "OC(=O)C");
}

TEST_F(TestStandardisation, TestChargedImidazole)
{
  _smiles = "CN1C=C[N+](CC)=C1";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _m2 = _m1;
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  EXPECT_EQ(_m1.molecular_formula(), "C6N2H11");
  EXPECT_EQ(_m1.smiles(), "CN1C=C[N+](=C1)CC");

  _chemical_standardisation.Activate(CS_CHARGED_IMIDAZOLE, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);

  _smiles = "CCN1C=C[N+](C)=C1";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_m1.molecular_formula(), "C6N2H11");
  EXPECT_EQ(_m1.unique_smiles(), "C[n+]1c[n](CC)cc1");

  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.unique_smiles(), "C[n]1c[n+](CC)cc1");

  _m1.invalidate_smiles();

  // The transformed molecule should not change

  constexpr int replicates = 10;
  for (int i = 0; i < replicates; ++i) {
    const IWString & smiles = _m1.random_smiles();
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smiles));
    _chemical_standardisation.process(m);
    EXPECT_EQ(m.unique_smiles(), "C[n]1c[n+](CC)cc1");
  }

  // Random variants should all get transformed to the same form.

  for (int i = 0; i < replicates; ++i) {
    const IWString & smiles = _m2.random_smiles();
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smiles));
    _chemical_standardisation.process(m);
    EXPECT_EQ(m.unique_smiles(), "C[n]1c[n+](CC)cc1");
  }
}

TEST_F(TestStandardisation, TestMisdrawnSulfonamideNoChange)
{
  _smiles = "CCS(=O)(=O)NC";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_MSDSA, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestMisdrawnSulfonamideChanges)
{
  _smiles = "CCS(=O)(O)=NC";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_MSDSA, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.unique_smiles(), "CNS(=O)(=O)CC");
}

TEST_F(TestStandardisation, TestEnoltoKetoYes)
{
  _smiles = "CC(O)=C";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_KETO_ENOL, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.unique_smiles(), "O=C(C)C");
}

TEST_F(TestStandardisation, TestEnoltKetoNoBecauseRing)
{
  _smiles = "C1C(O)=CC1";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_KETO_ENOL, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestEnoltKetoNoBecauseAdjacentKeto)
{
  _smiles = "CC(O)=CC(=O)C";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_KETO_ENOL, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestEnoltKetoNoBecauseAdjacentUnsaturation)
{
  _smiles = "C(=O)(O)C(=CC(=O)C1=CC=CC=C1OC)O";   //  CHEMBL4171775
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_CHARGED_IMIDAZOLE, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestEnoltKetoNoBecauseComplexInterdependency)
{
  _smiles = "C1=CC(=CC(=C1O)C(=O)CC(=O)C=CC1=CC=C(O)C=C1)Cl"; // CHEMBL4208282
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_CHARGED_IMIDAZOLE, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestEnoltKetoNoBecauseAdjacentUnsaturationDoubleBond)
{
  _smiles = "S(C1=NC2=CC(=CC=C2N1)C)CC(O)=C(C(=N)C)C#N"; // CHEMBL3197234
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_CHARGED_IMIDAZOLE, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestKetoToEnolNoTerminalDoubleBonds) {
  _smiles = "CC(=O)C";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _chemical_standardisation.set_keto_enol_preference(standardisation::KetoEnol::kToEnol);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);

  constexpr int kVerbose = false;

  _chemical_standardisation.Activate(CS_KETO_ENOL, kVerbose);

  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestKetoToEnolSimple) {
  _smiles = "CC(=O)CC";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _chemical_standardisation.set_keto_enol_preference(standardisation::KetoEnol::kToEnol);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);

  constexpr int kVerbose = false;

  _chemical_standardisation.Activate(CS_KETO_ENOL, kVerbose);

  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.unique_smiles(), "OC(=CC)C") << "from " << _smiles << " got "
                << _m1.unique_smiles();
}

TEST_F(TestStandardisation, TestKetoToEnolAdjacentNotProcessed) {
  _smiles = "OC(=O)CC(=O)CC(=O)C=CC(O)=O CHEMBL1743220";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _chemical_standardisation.set_keto_enol_preference(standardisation::KetoEnol::kToEnol);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);

  constexpr int kVerbose = false;

  _chemical_standardisation.Activate(CS_KETO_ENOL, kVerbose);

  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestKetoToEnolNotResolved) {
  _smiles = "CCC(=O)CC";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _chemical_standardisation.set_keto_enol_preference(standardisation::KetoEnol::kToEnol);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);

  constexpr int kVerbose = false;

  _chemical_standardisation.Activate(CS_KETO_ENOL, kVerbose);

  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestKetoToEnolHigherConnectivity) {
  _smiles = "CC(C)C(=O)CC";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _chemical_standardisation.set_keto_enol_preference(standardisation::KetoEnol::kToEnol);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);

  constexpr int kVerbose = false;

  _chemical_standardisation.Activate(CS_KETO_ENOL, kVerbose);

  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.unique_smiles(), "OC(=C(C)C)CC") << "from " << _smiles << " get " << 
            _m1.unique_smiles();
}

TEST_F(TestStandardisation, TestChargedPyrazole)
{
  _smiles = "[N+]1(=C(C)C=CN1CC1OC(=O)C(C1)(C1=CC=CC=C1)C1=CC=CC=C1)CC";  // CHEMBL140300
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _m2 = _m1;
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_CHARGED_IMIDAZOLE, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);

  EXPECT_EQ(_m1.unique_smiles(), "O=C1OC(CC1(c1ccccc1)c1ccccc1)C[n+]1[n](c(cc1)C)CC");

  _m1.invalidate_smiles();

  // The transformed molecule should not change

  constexpr int replicates = 10;
  for (int i = 0; i < replicates; ++i) {
    const IWString & smiles = _m1.random_smiles();
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smiles));
    _chemical_standardisation.process(m);
    EXPECT_EQ(m.unique_smiles(), "O=C1OC(CC1(c1ccccc1)c1ccccc1)C[n+]1[n](c(cc1)C)CC");
  }
  // Random variants of the starting molecule should all end up the same.
  for (int i = 0; i < replicates; ++i) {
    const IWString & smiles = _m2.random_smiles();
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smiles));
    _chemical_standardisation.process(m);
    EXPECT_EQ(m.unique_smiles(), "O=C1OC(CC1(c1ccccc1)c1ccccc1)C[n+]1[n](c(cc1)C)CC");
  }
}

// Since the pyrazole algorithm depends on the atom ordering, do the same test with
// different atom orderings.
TEST_F(TestStandardisation, TestChargedPyrazoleIncreasing) {
  _smiles = "ON1C=CC=[N+]1N";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_CHARGED_IMIDAZOLE, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);

  EXPECT_EQ(_m1.unique_smiles(), "O[n+]1[n](N)ccc1");
}

TEST_F(TestStandardisation, TestChargedPyrazoledecreasing) {
  _smiles = "N[N+]1=CC=CN1O";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_CHARGED_IMIDAZOLE, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);

  EXPECT_EQ(_m1.unique_smiles(), "O[n+]1[n](N)ccc1");
}

TEST_F(TestStandardisation, TestChargedImidazole3ConnectedNplus) {
  _smiles = "[O-][n]1cco[n+]1=C";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_CHARGED_IMIDAZOLE, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestHRemoval) {
  _smiles = "CC";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _m1.make_implicit_hydrogens_explicit();
  _chemical_standardisation.Activate(CS_XH, /*verbose*/ false);
  EXPECT_GT(_chemical_standardisation.process(_m1), 0);
}

#ifdef REDUNDANT_TEST
// Because of other transformations, this guard
// is never made.
TEST_F(TestStandardisation, NoChargedPyrazolones) {
  _smiles = "[O-]C1=C2C3=[N+](N1)C23";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _chemical_standardisation.activate_all();
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}
#endif

TEST_F(TestStandardisation, PyrazoloneNoProcess) {
  _smiles = "C1=N(=O)NC(=C1)O";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _chemical_standardisation.activate_all();
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

#ifdef THIS_CASE_DOES_NOT_EXIST
Jan 2024.
Looking at Chembl, this case seems not to exist. Turn off for now.
The existing standardisation does not change [N+]#[C-].
TEST_F(TestStandardisation, ReverseReversedCyano) {
  _smiles = "CCN#C";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  constexpr int kVerbose = 0;
  _chemical_standardisation.Activate(CS_REVERSE_NV5, kVerbose);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.unique_smiles(), "CC[N+]#[C-]") << _m1.unique_smiles() << " not match";
}
#endif

// First do a test without canonicalization, and then
// repeat with canonicalisation.
TEST_F(TestStandardisation, TestNoUsmi) {
  _smiles = "N1=C(N)C=CN=C1S";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _smiles = "N1=C(S)N=CC=C1N";
  ASSERT_TRUE(_m2.build_from_smiles(_smiles));

  _chemical_standardisation.activate_all();
  EXPECT_GT(_chemical_standardisation.process(_m1), 0);
  EXPECT_GT(_chemical_standardisation.process(_m2), 0);
  EXPECT_NE(_m1.unique_smiles(), _m2.unique_smiles());
}

// First do a test without canonicalization, and then
// repeat with canonicalisation.
TEST_F(TestStandardisation, TestUsmi) {
  _smiles = "N1=C(N)C=CN=C1S";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _smiles = "N1=C(S)N=CC=C1N";
  ASSERT_TRUE(_m2.build_from_smiles(_smiles));

  _chemical_standardisation.activate_all();
  _chemical_standardisation.set_convert_to_canonical_order(standardisation::Canonicalise::kReinterpretSmiles);
  EXPECT_GT(_chemical_standardisation.process(_m1), 0);
  EXPECT_GT(_chemical_standardisation.process(_m2), 0);
  EXPECT_EQ(_m1.unique_smiles(), _m2.unique_smiles());
}

TEST_F(TestStandardisation, TestCminus) {
  _smiles = "O=C1[C-](C(=O)C2=CC=CC=C12)C(=O)C1=CC=CC=C1";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _chemical_standardisation.activate_all();
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.unique_smiles(), "O=C(C1C(=O)c2c(C1=O)cccc2)c1ccccc1");
}

TEST_F(TestStandardisation, TestSulfonylUrea) {
  _smiles = "CNC(=N)S";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _chemical_standardisation.activate_all();
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.unique_smiles(), "S=C(NC)N") << 
        "sulfonyl urea not converted " << _m1.unique_smiles() << '\n';
}

TEST_F(TestStandardisation, TestSulfonylUreaRing) {
  _smiles = "N1=C(C(C)CN1C(S)=N)C1=CC=CC(=C1)Cl CHEMBL90428";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _chemical_standardisation.activate_all();
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.unique_smiles(), "Clc1cc(C2=NN(C(=S)N)CC2C)ccc1") << 
        "sulfonyl urea not converted " << _m1.unique_smiles() << '\n';
}

TEST_F(TestStandardisation, Test124Triazole) {
  _smiles = "SC1=NC(=C(C)N=N1)S CHEMBL3272977";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _chemical_standardisation.activate_all();
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.unique_smiles(), "Cc1[n][nH]c(=S)[nH]c1=S") << 
        "124 triazole not converted " << _m1.unique_smiles() << '\n';
}

TEST_F(TestStandardisation, TestExternalNoSmilesNoSmarts) {
  const_IWSubstring buffer = R"(
)";
  const_IWSubstring fname_not_used;

  EXPECT_FALSE(_chemical_standardisation.ReadExternalProto(buffer, fname_not_used));
}

TEST_F(TestStandardisation, TestExternalNoSmarts) {
  const_IWSubstring buffer = R"(
  smiles: "O=[N+]-[O-]"
)";
  const_IWSubstring fname_not_used;

  EXPECT_FALSE(_chemical_standardisation.ReadExternalProto(buffer, fname_not_used));
}

TEST_F(TestStandardisation, TestExternalNoSmiles) {
  const_IWSubstring buffer = R"(
  smarts: "O=N=O"
)";
  const_IWSubstring fname_not_used;

  EXPECT_FALSE(_chemical_standardisation.ReadExternalProto(buffer, fname_not_used));
}

TEST_F(TestStandardisation, TestExternalNitroToChargeSeparated) {
  const_IWSubstring buffer = R"(
  smiles: "O=[N+]-[O-]"
  smarts: "O=N=O"
)";
  const_IWSubstring fname_not_used;

  EXPECT_TRUE(_chemical_standardisation.ReadExternalProto(buffer, fname_not_used));
  ASSERT_TRUE(_chemical_standardisation.active());

  _smiles = "CN(=O)=O";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_TRUE(_chemical_standardisation.process(_m1));

  EXPECT_EQ(_m1.unique_smiles(), "O=[N+]([O-])C") << "GOt smiles " << _m1.smiles();
}

TEST_F(TestStandardisation, TestExternalChargeSeparatedToNitro) {
  const_IWSubstring buffer = R"(
  smarts: "O=[N+]-[O-]"
  smiles: "O=N=O"
)";
  const_IWSubstring fname_not_used;

  EXPECT_TRUE(_chemical_standardisation.ReadExternalProto(buffer, fname_not_used));
  ASSERT_TRUE(_chemical_standardisation.active());

  _smiles = "C[N+](=O)[O-]";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_TRUE(_chemical_standardisation.process(_m1));

  EXPECT_EQ(_m1.unique_smiles(), "CN(=O)=O") << "GOt smiles " << _m1.smiles();
}

TEST_F(TestStandardisation, TestExternalNeutraliseAcid) {
  const_IWSubstring buffer = R"(
  smarts: "[O-][C,S]=O"
  smiles: "O-*=O"
  name: "acid"
)";
  const_IWSubstring fname_not_used;

  EXPECT_TRUE(_chemical_standardisation.ReadExternalProto(buffer, fname_not_used));
  ASSERT_TRUE(_chemical_standardisation.active());
  _chemical_standardisation.set_append_string_depending_on_what_changed(1);

  _smiles = "[O-]C=O";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _m1.set_name("foo");
  EXPECT_TRUE(_chemical_standardisation.process(_m1));

  EXPECT_EQ(_m1.unique_smiles(), "OC=O") << "GOt smiles " << _m1.smiles();
  EXPECT_EQ(_m1.name(), "foo STD:acid") << "Name mismatch got " << _m1.name();
}

TEST_F(TestStandardisation, TestExternalChargedAcid) {
  const_IWSubstring buffer = R"(
  smarts: "[OH][C,S]=O"
  smiles: "[O-]-*=O"
  name: "acid"
)";
  const_IWSubstring fname_not_used;

  EXPECT_TRUE(_chemical_standardisation.ReadExternalProto(buffer, fname_not_used));
  ASSERT_TRUE(_chemical_standardisation.active());
  _chemical_standardisation.set_append_string_depending_on_what_changed(1);

  _smiles = "OC=O";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  _m1.set_name("foo");
  EXPECT_TRUE(_chemical_standardisation.process(_m1));

  EXPECT_EQ(_m1.unique_smiles(), "O=C[O-]") << "GOt smiles " << _m1.smiles();
  EXPECT_EQ(_m1.name(), "foo STD:acid") << "Name mismatch got " << _m1.name();
}

struct ExternalDirective {
  IWString directive;
  IWString smiles;
  int expected_rc;
  IWString expected_result;
};

class TestExternalDirective : public testing::TestWithParam<ExternalDirective> {
  protected:
    Chemical_Standardisation _chemical_standardisation;
    Molecule _m;
};

TEST_P(TestExternalDirective, Tests) {
  const auto params = GetParam();

  ASSERT_TRUE(_m.build_from_smiles(params.smiles));

  const_IWSubstring fname_not_used;
  EXPECT_TRUE(_chemical_standardisation.ReadExternalProto(params.directive, fname_not_used));
  EXPECT_EQ(_chemical_standardisation.process(_m), params.expected_rc);
  if (params.expected_rc > 0) {
    EXPECT_EQ(_m.unique_smiles(), params.expected_result) << " got " <<
                _m.unique_smiles() << " expected " << params.expected_result;
  }
}
INSTANTIATE_TEST_SUITE_P(TestExternalDirective, TestExternalDirective, testing::Values(
  ExternalDirective({R"(smarts: "[ND1H1]=[CD3]-[OD1H]" smiles: "N-C=O" name: "amide")",
                    "C1C2CCC1C(C2C(=N)O)C(=N)O CHEMBL1717199", 1,
                    "O=C(N)C1C(C(=O)N)C2CC1CC2"}),
  ExternalDirective({R"(smarts: "[NR0]=c1:[cD2]:[cD2]:[nH]:c:c1" smiles: "N-C=CC=N" name: "para-amino")",
                    "N(C(=NC(C)C)O)S(=O)(=O)C1=CNC=CC1=NC1=CC(C)=CC=C1 CHEMBL1148", 1,
                    "OC(=NC(C)C)NS(=O)(=O)c1c(Nc2cc(C)ccc2)cc[n]c1"})
));


struct ForStd {
  std::vector<IWString> directives;
  IWString smiles;
  IWString expected;
};


class TestStandardisationP: public testing::TestWithParam<ForStd> {
  protected:
    Chemical_Standardisation _chemical_standardisation;
    Molecule _m;
};

TEST_P(TestStandardisationP, Tests) {
  const auto params = GetParam();

  static constexpr int kVerbose = 0;

  for (const IWString& directive : params.directives) {
    ASSERT_TRUE(_chemical_standardisation.Activate(directive, kVerbose));
  }
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  ASSERT_TRUE(_chemical_standardisation.process(_m)) << "doing " << params.directives[0] <<
              ' ' << params.smiles;
  EXPECT_EQ(_m.unique_smiles(), params.expected) << "got " << 
            _m.unique_smiles() << " expected " << params.expected;

  for (int i = 0; i < 10; ++i) {
    _m.build_from_smiles(params.smiles);
    const IWString s = _m.random_smiles();
    ASSERT_TRUE(_m.build_from_smiles(s));
    ASSERT_TRUE(_chemical_standardisation.process(_m)) << "doing " << params.directives[0] <<
                ' ' << params.smiles << " rsmi " << s;
    EXPECT_EQ(_m.unique_smiles(), params.expected) << "got " << 
              _m.unique_smiles() << " expected " << params.expected << " rsmi " << s;
  }
}
INSTANTIATE_TEST_SUITE_P(TestStandardisationP, TestStandardisationP, testing::Values(
  ForStd{{"rvnv5"}, "N1(=NC(=N(=O)C2=CC(=CC=C12)OCCCN1CCOCC1)CC)=O CHEMBL553213", 
         "CCc1[n][n+]([O-])c2c([n+]1[O-])cc(OCCCN1CCOCC1)cc2"},
  ForStd{{"isotope"}, "[2H]-C", "C[H]"},
  ForStd{{"isotope", "all"}, "[2H]-C", "C"},

  ForStd{{"to2ap"}, "O=C(NN)c1cc[nH]c(=N)c1 CHEMBL3091876", "O=C(NN)c1cc([n]cc1)N"},
  ForStd{{"to2ap"}, "N=C(N)N=c1sc([n][nH]1)c1c(C)cccc1 CHEMBL1188079", "NC(=N)Nc1sc(c2c(C)cccc2)[n][n]1"},

  ForStd{{"isoxazole"}, "C12=C(SN=C1-O)CCCC2NC CHEMBL150489", "O=c1[nH]sc2c1C(NC)CCC2"},

  ForStd{{"oxopyrimidine"}, "C1(=C(NC(=NC1=O)C)NCC1=CN=CC=C1)C#N CHEMBL17125", "O=c1[nH]c(C)[n]c(NCc2c[n]ccc2)c1C#N"},
  ForStd{{"oxopyrimidine"}, "C1CCCC2=C1NC(=NC2=O)N=C(NC1=CC=C(C=C1)OC)N CHEMBL4932203",
                            "O=c1[nH]c([n]c2c1CCCC2)N=C(Nc1ccc(OC)cc1)N"},

//We are no longer standardising fused pyrazoles.
//ForStd{{"pyrazole"}, "C1=C2C(=NC(=NC2=NN1)C)O CHEMBL154781",
//              "Oc1c2c[n][nH]c2[n]c(C)[n]1"}

  ForStd{{"indoleh"}, "BrC1=CNC2=NC=NC2=C1 CHEMBL1562708",
                            "Brc1c[n]c2[n]c[nH]c2c1"},
  ForStd{{"indoleh"}, "CC1=NC(=C2C=CN=C2N1)NC(=O)[C@@H]1CCCOC1 CHEMBL5442493",
                      "O=C(Nc1c2c([nH]cc2)[n]c(C)[n]1)[C@H]1COCCC1"},
// Turned off because of the change to connectivity testing.
// Re-enable once that is resolved.
//ForStd{{"guan"}, "C(=NC1=CC=CN=C1)(NC#N)N1CCCCC1 CHEMBL86956",
//                          "N#CN=C(N1CCCCC1)Nc1c[n]ccc1"},

  ForStd{{"rvnv5o"}, "C1(=N(=O)C=CC=C1)O CHEMBL3898657",
                            "Oc1[n+]([O-])cccc1"},
  ForStd{{"rvnv5o"}, "N1(=C(C)CCC1C)=O CHEMBL332704",
                            "[O-][N+]1=C(CCC1C)C"}
));

struct ForStdKetoEnol {
  std::vector<IWString> directives;
  standardisation::KetoEnol keto_enol;
  IWString smiles;
  IWString expected;
};

class TestStandardisationKetoEnol : public testing::TestWithParam<ForStdKetoEnol> {
  protected:
    Chemical_Standardisation _chemical_standardisation;
    Molecule _m;
};

TEST_P(TestStandardisationKetoEnol, Tests) {
  const auto& params = GetParam();

  static constexpr int kVerbose = 1;

  for (const IWString& directive : params.directives) {
    ASSERT_TRUE(_chemical_standardisation.Activate(directive, kVerbose));
  }
  _chemical_standardisation.set_keto_enol_preference(params.keto_enol);

  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  EXPECT_TRUE(_chemical_standardisation.process(_m));

  EXPECT_EQ(_m.unique_smiles(), params.expected) << "got " << 
            _m.unique_smiles() << " expected " << params.expected << ' ' <<
            _m.name();

  ASSERT_TRUE(_m.valence_ok()) << "invalid valence " << _m.smiles() << ' ' <<
              _m.name();

  for (int i = 0; i < 10; ++i) {
    _m.build_from_smiles(params.smiles);
    const IWString s = _m.random_smiles();
    _m.build_from_smiles(s);
    ASSERT_TRUE(_chemical_standardisation.process(_m));
    EXPECT_EQ(_m.unique_smiles(), params.expected) << " from " << s << " got "
                << _m.unique_smiles() << " expected " << params.expected << '\n';
  }
}
INSTANTIATE_TEST_SUITE_P(TestStandardisationKetoEnol, TestStandardisationKetoEnol, testing::Values(
  ForStdKetoEnol{{"isoxazole"}, standardisation::KetoEnol::kToEnol, "C12=C(SNC1=O)CCCC2NC CHEMBL150489",
                 "Oc1[n]sc2c1C(NC)CCC2"},
  ForStdKetoEnol{{"isoxazole"}, standardisation::KetoEnol::kToKeto, "Oc1[n]sc2c1C(NC)CCC2 CHEMBL150489",
                 "O=c1[nH]sc2c1C(NC)CCC2"},

  ForStdKetoEnol{{"pyrazolone"}, standardisation::KetoEnol::kToEnol, "C1=CNNC1=O",
                 "Oc1[n][nH]cc1"},
  ForStdKetoEnol{{"pyrazolone"}, standardisation::KetoEnol::kToKeto, "Oc1[n][nH]cc1",
                 "O=c1[nH][nH]cc1"},
  ForStdKetoEnol{{"pyrazolone"}, standardisation::KetoEnol::kToEnol, "C1(=C(O)C(=O)NN1CCC)O CHEMBL366101.toEnol",
                "CCC[n]1[n]c(c(c1O)O)O"},
  ForStdKetoEnol{{"pyrazolone"}, standardisation::KetoEnol::kToKeto, "C1=CC(=CC=C1C1=NN(C)C(=C1)O)Cl CHEMBL4743178",
                "Clc1ccc(c2[nH][n](c(=O)c2)C)cc1"},
  ForStdKetoEnol{{"pyrazolone"}, standardisation::KetoEnol::kToKeto, "C(=O)(O)C1=CC=NC(=C1)N1N=CC=C1O CHEMBL4287317",
                "OC(=O)c1cc([n]2[nH]ccc2=O)[n]cc1"},
  ForStdKetoEnol{{"pyrazolone"}, standardisation::KetoEnol::kToKeto, "N1N=C(O)C=C1 CHEMBL4227850.toKeto",
                "O=c1[nH][nH]cc1"},
  ForStdKetoEnol{{"pyrazolone"}, standardisation::KetoEnol::kToEnol, "O=c1[nH][nH]cc1 CHEMBL4227850.ToEnol",
                "Oc1[n][nH]cc1"},
  ForStdKetoEnol{{"pyrazolone"}, standardisation::KetoEnol::kToKeto, "OC1=C(C)C(=NN1)C CHEMBL1476106.toKeto",
                "O=c1[nH][nH]c(c1C)C"},
  // Note that the round trip does not return to the starting molecule. We would then
  // need to run the molecule through pyrazole standardisation.
  ForStdKetoEnol{{"pyrazolone"}, standardisation::KetoEnol::kToEnol, "O=c1[nH][nH]c(c1C)C CHEMBL1476106.toEnol",
                "Oc1[n][nH]c(c1C)C"},

  ForStdKetoEnol{{"pyrazolone"}, standardisation::KetoEnol::kToKeto, "C1(=CC(=NN1)O)N toKeto",
                "O=c1cc(N)[nH][nH]1"},
  // Another transformation changes the exocyclic bond
  ForStdKetoEnol{{"pyrazolone", "enol-fused"}, standardisation::KetoEnol::kToKeto, "OC1=NN(C)C2=NC(C)=CC(C)=C12 toKeto",
                "O=c1[nH][n](c2[n]c(cc(c12)C)C)C"}


));

struct NoChange {
  std::vector<IWString> directives;
  IWString smiles;
};

class TestStandardisationNoChange : public testing::TestWithParam<NoChange> {
  protected:
    Chemical_Standardisation _chemical_standardisation;
    Molecule _m;
};

TEST_P(TestStandardisationNoChange, Tests) {
  const auto params = GetParam();

  static constexpr int kVerbose = 0;

  for (const IWString& directive : params.directives) {
    ASSERT_TRUE(_chemical_standardisation.Activate(directive, kVerbose));
  }
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  ASSERT_FALSE(_chemical_standardisation.process(_m)) << params.directives[0] << " changed " <<
                params.smiles;
}
INSTANTIATE_TEST_SUITE_P(TestStandardisationNoChange, TestStandardisationNoChange, testing::Values(
  NoChange{{"to2ap"}, "C12=NC(=N)C=CN1[C@H]1O[C@H]([C@H](O)[C@H]1O2)CO CHEMBL4303543"},
  NoChange{{"to2ap"}, "CC=CC1=C(C(=CC(=C1)CC1=CNC(=N)NC1=N)CCC)OC CHEMBL528943"},
  NoChange{{"indoleh"}, "CC1=CC(=O)N2N=CC=C2N1 CHEMBL4525057"},
  NoChange{{"guan"}, "C(=N)(NC1=C2C(=CC=C1)C=CC=C2)N1CCCCC1 CHEMBL91145"}
));

TEST_F(TestStandardisation, OxygenAndSulphurTautomersWithCombinedTransformations) {
  _chemical_standardisation.activate_all();
  // 4-pyridone is opt-in and is deliberately not enabled by activate_all().
  ASSERT_TRUE(_chemical_standardisation.Activate(CS_4_PYRIDONE, 0));
  for (const auto& [input, output] : std::vector<std::pair<const char*, const char*>>{
           {"CC([SH])=NC", "CC(=S)NC"},
           {"CC([OH])=[NH]", "CC(=O)N"},
           {"CC([SH])=[NH]", "CC(=S)N"},
           {"NC(S)=N", "NC(=S)N"},
           {"SC1=NCCCC1", "S=C1NCCCC1"},
           {"Sc1ccccn1", "S=C1C=CC=CN1"},
           {"CC([OH])=CC", "CC(=O)CC"},
           {"CC([SH])=CC", "CC(=S)CC"},
           {"[OH]c1cc[n]cc1", "O=C1C=CNC=C1"},
           {"[SH]c1cc[n]cc1", "S=C1C=CNC=C1"},
           {"O=c1cc[nH]c(C)n1", "O=c1ccnc(C)[nH]1"},
           {"S=c1cc[nH]c(C)n1", "S=c1ccnc(C)[nH]1"}}) {
    Molecule molecule;
    ASSERT_TRUE(molecule.build_from_smiles(input));
    Molecule expected;
    ASSERT_TRUE(expected.build_from_smiles(output));
    const IWString formula = molecule.molecular_formula();
    EXPECT_NE(_chemical_standardisation.process(molecule), 0);
    EXPECT_EQ(molecule.unique_smiles(), expected.unique_smiles()) << input;
    EXPECT_EQ(molecule.molecular_formula(), formula);
    EXPECT_TRUE(molecule.valence_ok());
    EXPECT_EQ(_chemical_standardisation.process(molecule), 0) << input;
  }
}

TEST_F(TestStandardisation, FourPyridoneRemainsOptIn) {
  _chemical_standardisation.activate_all();
  for (const char* input : {"Oc1ccncc1", "Sc1ccncc1"}) {
    Molecule molecule;
    ASSERT_TRUE(molecule.build_from_smiles(input));
    const IWString before = molecule.unique_smiles();
    EXPECT_EQ(_chemical_standardisation.process(molecule), 0);
    EXPECT_EQ(molecule.unique_smiles(), before);
  }
}

// Compare oxygen and sulfur tautomer rules, including formula preservation,
// valid valence, atom-order independence and idempotence.
struct TautomerCase {
  const char* directive;
  const char* input;
  const char* expected;
  bool changes;
  standardisation::KetoEnol preference = standardisation::KetoEnol::kToKeto;
};

class TestOxygenAndSulphurTautomers : public TestStandardisation,
                                 public testing::WithParamInterface<TautomerCase> {};

TEST_P(TestOxygenAndSulphurTautomers, CanonicalFormAndAtomOrder) {
  const TautomerCase& params = GetParam();
  ASSERT_TRUE(_chemical_standardisation.Activate(params.directive, 0));
  _chemical_standardisation.set_keto_enol_preference(params.preference);
  Molecule original;
  ASSERT_TRUE(original.build_from_smiles(params.input));
  Molecule expected;
  ASSERT_TRUE(expected.build_from_smiles(params.expected));
  const IWString expected_smiles = expected.unique_smiles();
  const IWString formula = original.molecular_formula();

  for (int i = 0; i < 10; ++i) {
    Molecule molecule;
    ASSERT_TRUE(molecule.build_from_smiles(i == 0 ? IWString(params.input) : original.random_smiles()));
    EXPECT_EQ(_chemical_standardisation.process(molecule) != 0, params.changes)
        << params.input;
    EXPECT_EQ(molecule.unique_smiles(), expected_smiles) << params.input;
    EXPECT_EQ(molecule.molecular_formula(), formula) << params.input;
    EXPECT_TRUE(molecule.valence_ok()) << molecule.smiles();
    EXPECT_EQ(_chemical_standardisation.process(molecule), 0)
        << "Standardisation must be idempotent: " << params.input;
  }
}

INSTANTIATE_TEST_SUITE_P(OxygenAndSulphur, TestOxygenAndSulphurTautomers, testing::Values(
    TautomerCase{CS_AMIDINE, "CC(N)=NC", "CC(=N)NC", true},
    TautomerCase{CS_AMIDINE, "CC([NH2])=[N]C", "CC(=N)NC", true},
    TautomerCase{CS_AMIDINE, "CC(=N)NC", "CC(=N)NC", false},
    TautomerCase{CS_AMIDINE, "C(N)=NC", "C(=N)NC", true},
    TautomerCase{CS_AMIDINE, "N=C(N)C", "N=C(N)C", false},
    TautomerCase{CS_AMIDINE, "CC(N)=NC1CCCCC1", "CC(=N)NC1CCCCC1", true},
    TautomerCase{CS_AMIDINE, "CC(N)=Nc1ccccc1", "CC(N)=Nc1ccccc1", false},
    TautomerCase{CS_AMIDINE, "CC(N)=NC(=N)N", "CC(N)=NC(=N)N", false},
    TautomerCase{CS_AMIDINE, "CC(N)=NNC", "CC(N)=NNC", false},
    TautomerCase{CS_AMIDINE, "CC(N)=NOC", "CC(N)=NOC", false},
    TautomerCase{CS_AMIDINE, "CC(N)=[NH+]C", "CC(N)=[NH+]C", false},
    TautomerCase{CS_AMIDINE, "C/C(N)=N/C", "C/C(N)=N/C", false},
    TautomerCase{CS_PYRIMIDINE, "N1=C(O)C=C(O)N=C1O", "N1C(=O)CC(=O)NC1=O", true},
    TautomerCase{CS_PYRIMIDINE, "N1C(=O)CC(=O)NC1=O", "N1C(=O)CC(=O)NC1=O", false},
    TautomerCase{CS_PYRIMIDINE, "N1C(=O)C=C(O)N=C1O", "N1C(=O)CC(=O)NC1=O", true},
    TautomerCase{CS_PYRIMIDINE, "N1=C(S)C=C(O)N=C1O", "N1C(=S)CC(=O)NC1=O", true},
    TautomerCase{CS_PYRIMIDINE, "N1=C(S)C=C(S)N=C1S", "N1C(=S)CC(=S)NC1=S", true},
    TautomerCase{CS_PYRIMIDINE, "CN1C(O)=CC(O)=NC1=O", "CN1C(=O)CC(=O)NC1=O", true},
    TautomerCase{CS_PYRIMIDINE, "CN1C(O)=CC(=O)N(C)C1=O", "CN1C(=O)CC(=O)N(C)C1=O", true},
    TautomerCase{CS_PYRIMIDINE, "N1=C(O)C(C)=C(O)N=C1O", "N1C(=O)C(C)C(=O)NC1=O", true},
    TautomerCase{CS_PYRIMIDINE, "N1=C(O)C(C)(C)C(O)=NC1=O", "N1C(=O)C(C)(C)C(=O)NC1=O", true},
    // Other ring arrangements, missing carbonyls, explicit isotope H and spiro
    // or marked chiral centers remain outside the coordinated family rule.
    TautomerCase{CS_PYRIMIDINE, "Oc1cc(O)ncn1", "Oc1cc(O)ncn1", false},
    TautomerCase{CS_PYRIMIDINE, "N1C(O)CC(O)NC1O", "N1C(O)CC(O)NC1O", false},
    TautomerCase{CS_PYRIMIDINE, "[2H]N1C(O)=CC(O)=NC1=O", "[2H]N1C(O)=CC(O)=NC1=O", false},
    TautomerCase{CS_PYRIMIDINE, "[NH2+]1C(=O)CC(=O)NC1=O", "[NH2+]1C(=O)CC(=O)NC1=O", false},
    TautomerCase{CS_PYRIMIDINE, "Oc1ncc(O)nc1O", "Oc1ncc(O)nc1O", false},
    TautomerCase{CS_PYRIMIDINE, "N1=C(O)C2(CC2)C(O)=NC1=O", "N1=C(O)C2(CC2)C(O)=NC1=O", false},
    TautomerCase{CS_PYRIMIDINE, "N1=C(O)[C@H](C)C(S)=NC1=O", "N1=C(O)[C@H](C)C(S)=NC1=O", false},
    TautomerCase{CS_PYRIMIDINE, "N1=C(O)C=C(O)N=C1O", "N1=C(O)C=C(O)N=C1O", false,
                 standardisation::KetoEnol::kToEnol},
    TautomerCase{CS_ALL, "CC(N)=NC", "CC(=N)NC", true},
    TautomerCase{CS_ALL, "N1=C(O)C=C(O)N=C1O", "N1C(=O)CC(=O)NC1=O", true},
    TautomerCase{CS_ALL, "N1=C(O)C=C(O)N=C1O.O=c1cc[nH]c(C)n1",
                 "N1C(=O)CC(=O)NC1=O.O=c1ccnc(C)[nH]1", true},
    TautomerCase{CS_ALL, "N1C(=O)CC(=O)NC1=O.O=c1cc[nH]c(C)n1",
                 "N1C(=O)CC(=O)NC1=O.O=c1ccnc(C)[nH]1", true},
    TautomerCase{CS_LACTIM_LACTAM, "CC(O)=NC", "CC(=O)NC", true},
    TautomerCase{CS_LACTIM_LACTAM, "CC(S)=NC", "CC(=S)NC", true},
    TautomerCase{CS_LACTIM_LACTAM, "CC([OH])=NC", "CC(=O)NC", true},
    TautomerCase{CS_LACTIM_LACTAM, "CC([SH])=NC", "CC(=S)NC", true},
    TautomerCase{CS_LACTIM_LACTAM, "CC(=O)NC", "CC(=O)NC", false},
    TautomerCase{CS_LACTIM_LACTAM, "CC(=S)NC", "CC(=S)NC", false},
    TautomerCase{CS_LACTIM_LACTAM_RING, "OC1=NCCCC1", "O=C1NCCCC1", true},
    TautomerCase{CS_LACTIM_LACTAM_RING, "SC1=NCCCC1", "S=C1NCCCC1", true},
    TautomerCase{CS_LACTIM_LACTAM_RING, "Oc1ccccn1", "O=C1C=CC=CN1", true},
    TautomerCase{CS_LACTIM_LACTAM_RING, "Sc1ccccn1", "S=C1C=CC=CN1", true},
    // Terminal neutral imine nitrogens are also supported.
    TautomerCase{CS_LACTIM_LACTAM, "CC(O)=N", "CC(=O)N", true},
    TautomerCase{CS_LACTIM_LACTAM, "CC(S)=N", "CC(=S)N", true},
    TautomerCase{CS_LACTIM_LACTAM, "CC([OH])=[NH]", "CC(=O)N", true},
    TautomerCase{CS_LACTIM_LACTAM, "CC([SH])=[NH]", "CC(=S)N", true},
    TautomerCase{CS_LACTIM_LACTAM, "NC(O)=N", "NC(=O)N", true},
    TautomerCase{CS_LACTIM_LACTAM, "NC(S)=N", "NC(=S)N", true},
    TautomerCase{CS_LACTIM_LACTAM, "CNC(O)=N", "CNC(=O)N", true},
    TautomerCase{CS_LACTIM_LACTAM, "CNC(S)=N", "CNC(=S)N", true},
    // Charged terminal imines and specified C=N stereochemistry stay excluded.
    TautomerCase{CS_LACTIM_LACTAM, "CC(O)=[NH2+]", "CC(O)=[NH2+]", false},
    TautomerCase{CS_LACTIM_LACTAM, "CC(S)=[NH2+]", "CC(S)=[NH2+]", false},
    TautomerCase{CS_LACTIM_LACTAM, "C/C(O)=N/C", "C/C(O)=N/C", false},
    TautomerCase{CS_LACTIM_LACTAM, "C/C(S)=N/C", "C/C(S)=N/C", false},
    // Two-connected central carbons and ether/thioether groups stay excluded.
    TautomerCase{CS_LACTIM_LACTAM, "OC=NC", "OC=NC", false},
    TautomerCase{CS_LACTIM_LACTAM, "SC=NC", "SC=NC", false},
    TautomerCase{CS_LACTIM_LACTAM, "COC(C)=NC", "COC(C)=NC", false},
    TautomerCase{CS_LACTIM_LACTAM, "CSC(C)=NC", "CSC(C)=NC", false},
    TautomerCase{CS_KETO_ENOL, "CC(O)=CC", "CC(=O)CC", true},
    TautomerCase{CS_KETO_ENOL, "CC(S)=CC", "CC(=S)CC", true},
    TautomerCase{CS_KETO_ENOL, "CC([OH])=CC", "CC(=O)CC", true},
    TautomerCase{CS_KETO_ENOL, "CC([SH])=CC", "CC(=S)CC", true},
    TautomerCase{CS_KETO_ENOL, "CC(=O)CC", "CC(=O)CC", false},
    TautomerCase{CS_KETO_ENOL, "CC(=S)CC", "CC(=S)CC", false},
    TautomerCase{CS_KETO_ENOL, "CC(=O)CC", "CC(O)=CC", true,
                 standardisation::KetoEnol::kToEnol},
    TautomerCase{CS_KETO_ENOL, "CC(=S)CC", "CC(S)=CC", true,
                 standardisation::KetoEnol::kToEnol},
    TautomerCase{CS_KETO_ENOL, "CC(=O)[CH2]C", "CC(O)=CC", true,
                 standardisation::KetoEnol::kToEnol},
    TautomerCase{CS_KETO_ENOL, "CC(=S)[CH2]C", "CC(S)=CC", true,
                 standardisation::KetoEnol::kToEnol},
    TautomerCase{CS_KETO_ENOL, "CC(=O)CC.CC(=O)CC", "CC(O)=CC.CC(O)=CC", true,
                 standardisation::KetoEnol::kToEnol},
    TautomerCase{CS_KETO_ENOL, "CC(=O)CC.CC(=S)CC", "CC(O)=CC.CC(S)=CC", true,
                 standardisation::KetoEnol::kToEnol},
    TautomerCase{CS_KETO_ENOL, "CC(=O)CCCCC(=S)C", "CC(O)=CCCC=C(S)C", true,
                 standardisation::KetoEnol::kToEnol},
    // Interacting keto/enol groups, including mixed O/S cases, remain excluded.
    TautomerCase{CS_KETO_ENOL, "CC(=O)C=C(S)C", "CC(=O)C=C(S)C", false},
    TautomerCase{CS_KETO_ENOL, "CC(=S)C=C(O)C", "CC(=S)C=C(O)C", false},
    TautomerCase{CS_KETO_ENOL, "CC(S)=CN", "CC(S)=CN", false},
    TautomerCase{CS_KETO_ENOL, "CC(SC)=CC", "CC(SC)=CC", false},
    TautomerCase{CS_KETO_ENOL, "OC1=CCCCC1", "OC1=CCCCC1", false},
    TautomerCase{CS_KETO_ENOL, "SC1=CCCCC1", "SC1=CCCCC1", false},
    TautomerCase{CS_KETO_ENOL, "CC(O)=CC.CC(S)=CC", "CC(=O)CC.CC(=S)CC", true},
    TautomerCase{CS_4_PYRIDONE, "Oc1ccncc1", "O=C1C=CNC=C1", true},
    TautomerCase{CS_4_PYRIDONE, "Sc1ccncc1", "S=C1C=CNC=C1", true},
    TautomerCase{CS_4_PYRIDONE, "[OH]c1cc[n]cc1", "O=C1C=CNC=C1", true},
    TautomerCase{CS_4_PYRIDONE, "[SH]c1cc[n]cc1", "S=C1C=CNC=C1", true},
    TautomerCase{CS_4_PYRIDONE, "O=C1C=CNC=C1", "O=C1C=CNC=C1", false},
    TautomerCase{CS_4_PYRIDONE, "S=C1C=CNC=C1", "S=C1C=CNC=C1", false},
    TautomerCase{CS_4_PYRIDONE, "CSc1ccncc1", "CSc1ccncc1", false},
    TautomerCase{CS_4_PYRIDONE, "Sc1cccnc1", "Sc1cccnc1", false},
    TautomerCase{CS_OXO_PYRIMIDINE, "O=c1cc[nH]c(C)n1", "O=c1ccnc(C)[nH]1", true},
    TautomerCase{CS_OXO_PYRIMIDINE, "S=c1cc[nH]c(C)n1", "S=c1ccnc(C)[nH]1", true},
    TautomerCase{CS_OXO_PYRIMIDINE, "O=c1ccnc(C)[nH]1", "O=c1ccnc(C)[nH]1", false},
    TautomerCase{CS_OXO_PYRIMIDINE, "S=c1ccnc(C)[nH]1", "S=c1ccnc(C)[nH]1", false},
    // Neither carbonyl nor thiocarbonyl competition may be ignored.
    TautomerCase{CS_OXO_PYRIMIDINE, "O=C1C2=C(C(=O)N=C(N1)N)N(C=N2)C1OC(CO)C(C1)O", "O=C1C2=C(C(=O)N=C(N1)N)N(C=N2)C1OC(CO)C(C1)O", false},
    TautomerCase{CS_OXO_PYRIMIDINE, "O=C1C2=C(C(=S)N=C(N1)N)N(C=N2)C1OC(CO)C(C1)O", "O=C1C2=C(C(=S)N=C(N1)N)N(C=N2)C1OC(CO)C(C1)O", false},
    TautomerCase{CS_OXO_PYRIMIDINE, "S=C1C2=C(C(=O)N=C(N1)N)N(C=N2)C1OC(CO)C(C1)O", "S=C1C2=C(C(=O)N=C(N1)N)N(C=N2)C1OC(CO)C(C1)O", false},
    TautomerCase{CS_OXO_PYRIMIDINE, "S=C1C2=C(C(=S)N=C(N1)N)N(C=N2)C1OC(CO)C(C1)O", "S=C1C2=C(C(=S)N=C(N1)N)N(C=N2)C1OC(CO)C(C1)O", false}
));

TEST_F(TestStandardisation, ProtonatedAmidineIsStableAfterOnePass) {
  _chemical_standardisation.activate_all();
  Molecule expected;
  ASSERT_TRUE(expected.build_from_smiles("CC(=N)NC"));
  for (int i = 0; i < 100; ++i) {
    Molecule initial;
    ASSERT_TRUE(initial.build_from_smiles("CC(N)=[NH+]C"));
    Molecule molecule;
    ASSERT_TRUE(molecule.build_from_smiles(initial.random_smiles()));
    EXPECT_GT(_chemical_standardisation.process(molecule), 0);
    EXPECT_EQ(molecule.unique_smiles(), expected.unique_smiles());
    EXPECT_TRUE(molecule.valence_ok());
    EXPECT_EQ(_chemical_standardisation.process(molecule), 0);
  }
}

TEST_F(TestStandardisation, NewTautomerRulesCanBeDisabled) {
  _chemical_standardisation.activate_all();
  ASSERT_TRUE(_chemical_standardisation.Activate("-amidine", 0));
  ASSERT_TRUE(_chemical_standardisation.Activate("-pyrimidine", 0));
  for (const char* input : {"CC(N)=NC", "N1=C(O)C=C(O)N=C1O"}) {
    Molecule molecule;
    ASSERT_TRUE(molecule.build_from_smiles(input));
    const IWString before = molecule.unique_smiles();
    EXPECT_EQ(_chemical_standardisation.process(molecule), 0);
    EXPECT_EQ(molecule.unique_smiles(), before);
  }
}

}  // namespace
