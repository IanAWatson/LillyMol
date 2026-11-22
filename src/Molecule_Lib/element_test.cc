#include <initializer_list>
#include <memory>
#include <limits>
#include <set>
#include <vector>

#include "gtest/gtest.h"

#include "Foundational/cmdline/cmdline.h"
#include "element.h"

namespace {

class ElementOptionsTest : public testing::Test {
 protected:
  int _selenium_organic;
  int _hydrogen_brackets;
  void SetUp() override {
    reset_element_file_scope_variables();
    _selenium_organic = get_element_from_atomic_number(34)->organic();
    _hydrogen_brackets = get_element_from_atomic_number(1)->needs_square_brackets();
  }
  void TearDown() override {
    const_cast<Element*>(get_element_from_atomic_number(34))->set_organic(_selenium_organic);
    set_explicit_hydrogens_need_square_brackets_in_smiles(_hydrogen_brackets);
    reset_element_file_scope_variables();
  }
  int Process(std::initializer_list<const char*> args) {
    std::vector<char*> argv{const_cast<char*>("element_test")};
    for (const char* arg : args) {
      argv.push_back(const_cast<char*>("-E"));
      argv.push_back(const_cast<char*>(arg));
    }
    Command_Line cl(argv.size(), argv.data(), "E:");
    return process_elements(cl);
  }
};

TEST_F(ElementOptionsTest, MarkExistingElementOrganic) {
  ASSERT_TRUE(Process({"O:Se"}));
  EXPECT_TRUE(get_element_from_atomic_number(34)->organic());
}

TEST_F(ElementOptionsTest, ExplicitShortNameNeedsNoAutocreate) {
  ASSERT_TRUE(Process({"Zz"}));
  const Element* element = get_element_from_symbol_no_case_conversion("Zz");
  ASSERT_NE(element, nullptr);
  EXPECT_TRUE(element->ok());
  EXPECT_FALSE(element->organic());
}

TEST_F(ElementOptionsTest, ModifierCanCreateShortName) {
  ASSERT_TRUE(Process({"O:Zz"}));
  const Element* element = get_element_from_symbol_no_case_conversion("Zz");
  ASSERT_NE(element, nullptr);
  EXPECT_TRUE(element->ok());
  EXPECT_TRUE(element->organic());
}

TEST_F(ElementOptionsTest, LongOrganicNameRequiresAnyLength) {
  EXPECT_FALSE(Process({"O:Foo"}));
  EXPECT_EQ(get_element_from_symbol_no_case_conversion("Foo"), nullptr);
}

TEST_F(ElementOptionsTest, AutocreateAloneDoesNotEnableLongModifiers) {
  EXPECT_FALSE(Process({"autocreate", "O:Foo"}));
  EXPECT_EQ(get_element_from_symbol_no_case_conversion("Foo"), nullptr);
}

TEST_F(ElementOptionsTest, LongOrganicNameWithAnyLength) {
  ASSERT_TRUE(Process({"anylength", "O:Foo"}));
  const Element* element = get_element_from_symbol_no_case_conversion("Foo");
  ASSERT_NE(element, nullptr);
  EXPECT_TRUE(element->ok());
  EXPECT_TRUE(element->organic());
  EXPECT_EQ(element->normal_valence(), kValenceNotDefined);
}

TEST_F(ElementOptionsTest, ExplicitLongNameDefaultsToNonOrganic) {
  ASSERT_TRUE(Process({"anylength", "Foo"}));
  const Element* element = get_element_from_symbol_no_case_conversion("Foo");
  ASSERT_NE(element, nullptr);
  EXPECT_TRUE(element->ok());
  EXPECT_FALSE(element->organic());
}

TEST_F(ElementOptionsTest, EmptyModifierNamesFail) {
  EXPECT_FALSE(Process({"O:"}));
  EXPECT_FALSE(Process({"nsqb="}));
}

TEST_F(ElementOptionsTest, LongBracketModifierRequiresAnyLength) {
  EXPECT_FALSE(Process({"nsqb=Foo"}));
  EXPECT_FALSE(Process({"autocreate", "nsqb=Foo"}));
  EXPECT_EQ(get_element_from_symbol_no_case_conversion("Foo"), nullptr);
}

TEST_F(ElementOptionsTest, LongBracketModifierWithAnyLength) {
  ASSERT_TRUE(Process({"anylength", "nsqb=Foo"}));
  const Element* element = get_element_from_symbol_no_case_conversion("Foo");
  ASSERT_NE(element, nullptr);
  EXPECT_TRUE(element->ok());
  EXPECT_TRUE(element->needs_square_brackets());
  EXPECT_FALSE(element->organic());
}

TEST_F(ElementOptionsTest, AnyLengthMustPrecedeLongName) {
  EXPECT_FALSE(Process({"O:Foo", "anylength"}));
  EXPECT_EQ(get_element_from_symbol_no_case_conversion("Foo"), nullptr);
}

TEST_F(ElementOptionsTest, AutocreateAloneCannotRegisterInvalidLongElement) {
  set_auto_create_new_elements(1);
  EXPECT_EQ(create_element_with_symbol("Foo"), nullptr);
  EXPECT_EQ(get_element_from_symbol_no_case_conversion("Foo"), nullptr);
}

TEST_F(ElementOptionsTest, ModifierCanUseExplicitlyDeclaredLongElement) {
  ASSERT_TRUE(Process({"anylength", "Foo", "O:Foo", "nsqb=Foo"}));
  const Element* element = get_element_from_symbol_no_case_conversion("Foo");
  ASSERT_NE(element, nullptr);
  EXPECT_TRUE(element->organic());
  EXPECT_TRUE(element->needs_square_brackets());
}
TEST_F(ElementOptionsTest, HydrogenBracketsCanBeDisabledAndRestored) {
  const Element* hydrogen = get_element_from_atomic_number(1);
  set_explicit_hydrogens_need_square_brackets_in_smiles(0);
  EXPECT_FALSE(hydrogen->needs_square_brackets());
  set_explicit_hydrogens_need_square_brackets_in_smiles(1);
  EXPECT_TRUE(hydrogen->needs_square_brackets());
}

TEST_F(ElementOptionsTest, IsotopeOnlyBufferDoesNotReadPastEnd) {
  // An exactly sized allocation makes an extra read detectable under ASan.
  auto buffer = std::make_unique<char[]>(2);
  buffer[0] = '1';
  buffer[1] = '2';
  isotope_t isotope;
  EXPECT_EQ(get_element_from_symbol(buffer.get(), 2, isotope), nullptr);
}

TEST_F(ElementOptionsTest, BoundedIsotopeAndElementNeedNoTerminator) {
  const char buffer[] = {'1', '3', 'C'};
  isotope_t isotope;
  EXPECT_EQ(get_element_from_symbol(buffer, 3, isotope), get_element_from_atomic_number(6));
  EXPECT_EQ(isotope, 13u);
}

TEST_F(ElementOptionsTest, MaximumIsotopeAndOverflow) {
  isotope_t isotope;
  EXPECT_EQ(get_element_from_symbol("4294967295C", isotope), get_element_from_atomic_number(6));
  EXPECT_EQ(isotope, std::numeric_limits<isotope_t>::max());
  EXPECT_EQ(get_element_from_symbol("4294967296C", isotope), nullptr);
  EXPECT_EQ(isotope, 0u);
  EXPECT_EQ(get_element_from_symbol("9999999999999999999999999999C", isotope), nullptr);
}

TEST_F(ElementOptionsTest, EmptyBoundedSymbolFails) {
  isotope_t isotope = 13;
  EXPECT_EQ(get_element_from_symbol(nullptr, 0, isotope), nullptr);
  EXPECT_EQ(isotope, 0u);
  EXPECT_EQ(get_element_from_symbol("C", -1, isotope), nullptr);
  EXPECT_EQ(get_element_from_symbol_no_case_conversion(nullptr, 0), nullptr);
}

TEST_F(ElementOptionsTest, InvalidSymbolsDoNotBecomeDummyAtoms) {
  isotope_t isotope;
  for (const char* symbol : {"!", "@", "C!", "?"}) {
    EXPECT_EQ(get_element_from_symbol_no_case_conversion(symbol), nullptr) << symbol;
    EXPECT_EQ(get_element_from_symbol(symbol, isotope), nullptr) << symbol;
  }
  EXPECT_EQ(get_element_from_symbol_no_case_conversion("*"), get_element_from_atomic_number(0));
}

TEST_F(ElementOptionsTest, InvalidSmartsSymbolWithAutocreateDoesNotIndexNegativeHash) {
  set_auto_create_new_elements(1);
  const Element* result = get_element_from_atomic_number(6);
  EXPECT_EQ(element_from_smarts_string("!", 1, result), 0);
  EXPECT_EQ(result, nullptr);
  EXPECT_EQ(element_from_smarts_string("*", 1, result), 1);
  EXPECT_EQ(result, get_element_from_atomic_number(0));
}

TEST_F(ElementOptionsTest, EmptySmilesAndSmartsSymbolsFail) {
  const Element* result = get_element_from_atomic_number(6);
  EXPECT_EQ(element_from_smiles_string(nullptr, 0, result), 0);
  EXPECT_EQ(result, nullptr);
  EXPECT_EQ(element_from_smarts_string(nullptr, 0, result), 0);
  EXPECT_EQ(result, nullptr);
  EXPECT_EQ(element_from_long_smiles_string(nullptr, 0, result), 0);
  EXPECT_EQ(result, nullptr);
}

TEST_F(ElementOptionsTest, OrdinarySmilesAndSmartsElementsStillWork) {
  const Element* result;
  EXPECT_EQ(element_from_smiles_string("Cl", 2, result), 2);
  EXPECT_EQ(result, get_element_from_atomic_number(17));
  EXPECT_EQ(element_from_smarts_string("CD3", 3, result), 1);
  EXPECT_EQ(result, get_element_from_atomic_number(6));
  EXPECT_EQ(element_from_smarts_string("nH", 2, result), 1);
  EXPECT_EQ(result, get_element_from_atomic_number(7));
}

TEST_F(ElementOptionsTest, PeriodicElementIdsRemainStable) {
  for (int z = 0; z <= HIGHEST_ATOMIC_NUMBER; ++z) {
    EXPECT_EQ(get_element_from_atomic_number(z)->unique_id(), z + 1);
  }
  Element carbon(6);
  EXPECT_EQ(carbon.unique_id(), get_element_from_atomic_number(6)->unique_id());
}

TEST_F(ElementOptionsTest, DirectAndRegisteredCustomElementsHaveDistinctIds) {
  Element direct("Zz");
  ASSERT_TRUE(Process({"Yy", "anylength", "LongCustomElement"}));
  std::set<int> ids;
  ids.insert(direct.unique_id());
  ids.insert(get_element_from_symbol_no_case_conversion("Yy")->unique_id());
  ids.insert(get_element_from_symbol_no_case_conversion("LongCustomElement")->unique_id());
  EXPECT_EQ(ids.size(), 3u);
  EXPECT_GT(*ids.begin(), HIGHEST_ATOMIC_NUMBER + 1);
}

TEST_F(ElementOptionsTest, CustomIdsAreNotReusedAcrossReset) {
  Element direct("Zz");
  ASSERT_TRUE(Process({"Yy"}));
  const int old_id = get_element_from_symbol_no_case_conversion("Yy")->unique_id();
  reset_element_file_scope_variables();
  EXPECT_EQ(get_element_from_symbol_no_case_conversion("Yy"), nullptr);
  EXPECT_EQ(get_element_from_symbol_no_case_conversion("Zz"), &direct);
  ASSERT_TRUE(Process({"Yy"}));
  const int new_id = get_element_from_symbol_no_case_conversion("Yy")->unique_id();
  EXPECT_NE(new_id, old_id);
  EXPECT_NE(new_id, direct.unique_id());
}

TEST_F(ElementOptionsTest, DirectElementDestructionClearsItsHashEntry) {
  {
    Element direct("Zz");
    EXPECT_EQ(get_element_from_symbol_no_case_conversion("Zz"), &direct);
  }
  EXPECT_EQ(get_element_from_symbol_no_case_conversion("Zz"), nullptr);
}

TEST_F(ElementOptionsTest, DestructionDoesNotClearAnotherElementsHashEntry) {
  const Element* carbon = get_element_from_atomic_number(6);
  {
    Element another_carbon(6);
    EXPECT_EQ(get_element_from_symbol_no_case_conversion("C"), carbon);
  }
  EXPECT_EQ(get_element_from_symbol_no_case_conversion("C"), carbon);
}

TEST_F(ElementOptionsTest, ResetClearsCustomHashesButPreservesPeriodicElements) {
  const Element* carbon = get_element_from_atomic_number(6);
  ASSERT_TRUE(Process({"Zz", "anylength", "LongCustomElement"}));
  reset_element_file_scope_variables();
  EXPECT_EQ(get_element_from_symbol_no_case_conversion("Zz"), nullptr);
  EXPECT_EQ(get_element_from_symbol_no_case_conversion("LongCustomElement"), nullptr);
  EXPECT_EQ(get_element_from_symbol_no_case_conversion("C"), carbon);
}

TEST_F(ElementOptionsTest, HashDecoderReplacesDestinationForEverySymbolType) {
  ASSERT_TRUE(Process({"??", "R#", "Zz", "R9", "a"}));
  for (const char* symbol : {"C", "Cl", "*", "??", "R#", "Zz", "R9", "a"}) {
    const Element* element = get_element_from_symbol_no_case_conversion(symbol);
    ASSERT_NE(element, nullptr) << symbol;
    IWString decoded("previous content");
    EXPECT_TRUE(symbol_for_atomic_symbol_hash_value(element->atomic_symbol_hash_value(), decoded));
    EXPECT_EQ(decoded, symbol);
  }
}

TEST_F(ElementOptionsTest, HashDecoderRejectsOutOfRangeValues) {
  IWString decoded("unchanged");
  EXPECT_FALSE(symbol_for_atomic_symbol_hash_value(-1, decoded));
  EXPECT_FALSE(symbol_for_atomic_symbol_hash_value(std::numeric_limits<int>::max(), decoded));
  EXPECT_EQ(decoded, "unchanged");
}

// Deallocation is terminal, so exercise it in a child process only.
TEST(ElementCleanupTest, DeallocationClearsHashLookup) {
  EXPECT_EXIT({
    de_allocate_periodic_table();
    const bool empty = get_element_from_symbol_no_case_conversion("C") == nullptr &&
                       get_element_from_symbol_no_case_conversion("*") == nullptr;
    std::exit(empty ? 0 : 1);
  }, testing::ExitedWithCode(0), "");
}

}  // namespace
