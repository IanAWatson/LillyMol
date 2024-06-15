// Parse smarts.

#include <iostream>

#include "Foundational/iwmisc/misc.h"

#include "misc2.h"
#include "molecule.h"
#include "parse_smarts_tmp.h"
#include "smiles.h"
#include "substructure.h"

using std::cerr;

using smiles::fetch_ring_number;

// Duplicated from parse_smiles. Just easier this way.

// For the Organic Subset, we need rapid access to certain elements.

static const Element * smi_element_star = nullptr;   // initialised for check in parse_smiles_token
static const Element * smi_element_b;
static const Element * smi_element_c;
static const Element * smi_element_n;
static const Element * smi_element_o;
static const Element * smi_element_f;
static const Element * smi_element_p;
static const Element * smi_element_s;
static const Element * smi_element_cl;
static const Element * smi_element_br;
static const Element * smi_element_i;
static const Element * smi_element_hydrogen;

static void
initialise_organic_subset()
{
  smi_element_star     = get_element_from_atomic_number(0);
  smi_element_hydrogen = get_element_from_atomic_number(1);
  smi_element_b    = get_element_from_atomic_number(5);
  smi_element_c    = get_element_from_atomic_number(6);
  smi_element_n    = get_element_from_atomic_number(7);
  smi_element_o    = get_element_from_atomic_number(8);
  smi_element_f    = get_element_from_atomic_number(9);
  smi_element_p    = get_element_from_atomic_number(15);
  smi_element_s    = get_element_from_atomic_number(16);
  smi_element_cl   = get_element_from_atomic_number(17);
  smi_element_br   = get_element_from_atomic_number(35);
  smi_element_i    = get_element_from_atomic_number(53);

  return;
}

// compatability table for tokens in a smarts.

#define PREVIOUS_TOKEN_WAS_NOT_SPECIFIED 0
#define PREVIOUS_TOKEN_WAS_RING 1
#define PREVIOUS_TOKEN_WAS_BOND 2
#define PREVIOUS_TOKEN_WAS_ATOM 3
#define PREVIOUS_TOKEN_WAS_OPEN_PAREN 4
#define PREVIOUS_TOKEN_WAS_CLOSE_PAREN 5
#define PREVIOUS_TOKEN_WAS_THREE_DOTS 6

#define CTDIM 7

static int * compatability_table = nullptr;     // never freed

static int
initialise_compatability_table()
{
  compatability_table = new_int(CTDIM * CTDIM);

// What can begin a smarts

  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_NOT_SPECIFIED + PREVIOUS_TOKEN_WAS_ATOM] = 1;

// What can follow a ring

  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_RING + PREVIOUS_TOKEN_WAS_RING] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_RING + PREVIOUS_TOKEN_WAS_BOND] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_RING + PREVIOUS_TOKEN_WAS_ATOM] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_RING + PREVIOUS_TOKEN_WAS_OPEN_PAREN] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_RING + PREVIOUS_TOKEN_WAS_CLOSE_PAREN] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_RING + PREVIOUS_TOKEN_WAS_THREE_DOTS] = 1;

// What can follow a bond

  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_BOND + PREVIOUS_TOKEN_WAS_RING] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_BOND + PREVIOUS_TOKEN_WAS_ATOM] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_BOND + PREVIOUS_TOKEN_WAS_OPEN_PAREN] = 1;

// What can follow an atom

  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_ATOM + PREVIOUS_TOKEN_WAS_ATOM] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_ATOM + PREVIOUS_TOKEN_WAS_BOND] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_ATOM + PREVIOUS_TOKEN_WAS_OPEN_PAREN] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_ATOM + PREVIOUS_TOKEN_WAS_CLOSE_PAREN] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_ATOM + PREVIOUS_TOKEN_WAS_RING] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_ATOM + PREVIOUS_TOKEN_WAS_THREE_DOTS] = 1;

// What can follow an open paren

  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_OPEN_PAREN + PREVIOUS_TOKEN_WAS_RING] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_OPEN_PAREN + PREVIOUS_TOKEN_WAS_BOND] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_OPEN_PAREN + PREVIOUS_TOKEN_WAS_ATOM] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_OPEN_PAREN + PREVIOUS_TOKEN_WAS_THREE_DOTS] = 1;

// What can follow a close paren

  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_CLOSE_PAREN + PREVIOUS_TOKEN_WAS_RING] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_CLOSE_PAREN + PREVIOUS_TOKEN_WAS_BOND] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_CLOSE_PAREN + PREVIOUS_TOKEN_WAS_ATOM] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_CLOSE_PAREN + PREVIOUS_TOKEN_WAS_OPEN_PAREN] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_CLOSE_PAREN + PREVIOUS_TOKEN_WAS_CLOSE_PAREN] = 1;
  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_CLOSE_PAREN + PREVIOUS_TOKEN_WAS_THREE_DOTS] = 1;

// What can follow three dots

  compatability_table[CTDIM * PREVIOUS_TOKEN_WAS_THREE_DOTS + PREVIOUS_TOKEN_WAS_ATOM] = 1;

  return 1;
}

static int
check_compatiability_table(int & previous_token_was, int nt)
{
  if (0 == compatability_table[CTDIM * previous_token_was + nt])
    return 0;

  previous_token_was = nt;

  return 1;
}

static int
is_three_dots(const const_IWSubstring & smarts,
              int characters_processed,
              const_IWSubstring & three_dots_qualifier)
{
//cerr << "Checking for three dots, in '" << smarts << "', characters_processed = " << characters_processed << " examine '" << smarts[characters_processed] << "'\n";
  if ('.' != smarts[characters_processed])
    return 0;

  int dots = 1;

  int smarts_length = smarts.length();

  characters_processed++;

  while (characters_processed < smarts_length)
  {
    if ('.' == smarts[characters_processed])
    {
      dots++;
      characters_processed++;
    }
    else
      break;
  }

  if (3 != dots)
    return 0;

  three_dots_qualifier.make_empty();

  if (characters_processed == smarts_length)   // wrong, smarts cannot end with ...
    return 1;

  if ('{' != smarts[characters_processed])    // no qualifier
    return 1;

  for (int i = characters_processed + 1; i < smarts_length; i++)
  {
    if ('}' != smarts[i])
      continue;

    smarts.from_to(characters_processed, i, three_dots_qualifier);
    return 1;
  }

// Never found a closing brace, bad news

  smarts.from_to(characters_processed + 3, smarts_length - 1, three_dots_qualifier);

  return 1;
}

/*
  In order to parse the '...' directive, we need to tell our caller what was
  the identity of the last atom created
*/

//#define DEBUG_BUILD_FROM_SMARTS

#define LOOKS_LIKE_SMARTS_BOND(b) ((SINGLE_BOND_SYMBOL == (b)) || (DOUBLE_BOND_SYMBOL == (b)) || (TRIPLE_BOND_SYMBOL == (b)) ||\
                                   (AROMATIC_BOND_SYMBOL == (b)) || ('~' == (b)) || ('@' == (b)) ||\
                                   ('!' == (b)) || ('&' == (b)) || (',' == (b)) ||\
                                   (';' == (b)) || ('^' == (b)) )

/*
  LAST_ATOM_CREATED keeps track of the atom numbers in the current fragment.
  ATOMS_IN_PREVIOUS_DISCONNECTED_SECTIONS is needed in order to be able to assign value
  initial atom numbers

*/

int
Substructure_Atom::_parse_smarts_specifier(const const_IWSubstring & qsmarts,
                                Parse_Smarts_Tmp & pst,
                                int atoms_in_previous_disconnected_sections,
                                Smiles_Ring_Status & ring_status)
{
#ifdef DEBUG_BUILD_FROM_SMARTS
  cerr << "Substructure_Atom::_parse_smarts_specifier: '" << qsmarts << "', " << atoms_in_previous_disconnected_sections << " atoms in previous fragments\n";
#endif

  extending_resizable_array<Substructure_Atom *> & completed = pst.completed();

  int characters_to_process = qsmarts.nchars();
  const char * smarts = qsmarts.rawchars();

  if (nullptr == smi_element_star) {     // initialise elements first time through
    initialise_organic_subset();
  }

  if (nullptr == compatability_table) {
    initialise_compatability_table();
  }

  int characters_processed = 0;

// We need to be somewhat careful about resizing, as we may be called recursively

// Various stacks to keep track of branches. Probably should be combined into
// a single kind with a new object.

  resizable_array<Substructure_Atom *>  atom_stack;
  resizable_array<int>                  chirality_stack;

  int previous_token_was = PREVIOUS_TOKEN_WAS_NOT_SPECIFIED;

  std::unique_ptr<Substructure_Bond> previous_bond;
//Substructure_Bond * previous_bond = nullptr;
  Substructure_Atom * previous_atom = nullptr;
  int previous_atom_chiral_count = 0;

  int atoms_this_fragment = 0;

  int paren_level = 0;

  const_IWSubstring three_dots_qualifier;    // any qualifier after a ... directive

  while (characters_processed < characters_to_process)
  {
    const char * s = smarts + characters_processed;

#ifdef DEBUG_BUILD_FROM_SMARTS
    cerr << "Examining smarts token '" << *s << "', previous " << previous_token_was << '\n';
#endif

    if (characters_processed && (isdigit(*s) || '%' == *s))    // we do not properly handle chirality here, fix sometime...
    {
      int ring_number;
      int nchars;
      if (! fetch_ring_number(s, characters_to_process - characters_processed, ring_number, nchars))
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Invalid ring specification");
        return 0;
      }

#ifdef DEBUG_BUILD_FROM_SMARTS
      cerr << "Is ring number " << ring_number << '\n';
#endif

      std::unique_ptr<Substructure_Bond> b;
      bond_type_t bt = SINGLE_BOND; 
      if (PREVIOUS_TOKEN_WAS_BOND == previous_token_was)
      {
        b.reset(previous_bond.release());
        bt = b->types_matched();
      }
      else
      {
        b.reset(new Substructure_Bond);
        b->make_single_or_aromatic();
        bt = NOT_A_BOND;
      }

      atom_number_t other_end;    // will be set when a ring closure
      if (ring_status.encounter(ring_number, previous_atom->unique_id(), other_end, bt))   // ring closing.
      {
        assert(nullptr != completed[other_end]);

#ifdef DEBUG_BUILD_FROM_SMARTS
        cerr << "Ring closure, atom at other end is " << other_end << " prev " << previous_atom->unique_id() << '\n';
#endif

        if (previous_atom->is_bonded_to(other_end))
        {
          cerr << "Substructure_Atom::_parse_smarts_specifier:two membered ring encountered\n";
          return 0;
        }

        b->set_atom(completed[other_end]);
        b->set_type(bt);
        previous_atom->_add_bond(b.release());
      }

      characters_processed += nchars;
      previous_token_was = PREVIOUS_TOKEN_WAS_RING;
    }
    else if (LOOKS_LIKE_SMARTS_BOND(*s) && check_compatiability_table(previous_token_was, PREVIOUS_TOKEN_WAS_BOND))
    {
      if (previous_bond)
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Consecutive bonds");
        return 0;
      }

      previous_bond.reset(new Substructure_Bond);
      int ncb;
      if (! previous_bond->construct_from_smarts(s, characters_to_process - characters_processed, ncb))
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Cannot parse bond specifier");
        return 0;
      }

      characters_processed += ncb;
#ifdef DEBUG_BUILD_FROM_SMARTS
      cerr << "Recognised as bond\n";
#endif
    }
    else if (is_three_dots(qsmarts, characters_processed, three_dots_qualifier) && check_compatiability_table(previous_token_was, PREVIOUS_TOKEN_WAS_THREE_DOTS))
    {
      if (previous_bond)
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Cannot have ... after bond");
        return 0;
      }

#ifdef DEBUG_BUILD_FROM_SMARTS
      cerr << "Got three dots, qualifier '" << three_dots_qualifier << "'\n";
#endif

      characters_processed += 3 + three_dots_qualifier.length();

      ThreeDots * three_dots = new ThreeDots(previous_atom->unique_id(),
                                             pst.last_query_atom_created() + 1);

      if (three_dots_qualifier.length())
        three_dots->set_qualifier(three_dots_qualifier);

      pst.add_three_dots(three_dots);

      Substructure_Atom * a = new Substructure_Atom;
      pst.add_root_atom(a);

      const_IWSubstring newsmarts(qsmarts);
      newsmarts.remove_leading_chars(characters_processed);

      if (0 == newsmarts.length())
      {
        cerr << "Substructure_Atom::parse_smiles_token:smarts cannot end in ... operator\n";
        delete a;
        return 0;
      }

      return a->_parse_smarts_specifier(newsmarts, pst, atoms_in_previous_disconnected_sections + atoms_this_fragment, ring_status);
    }
    else if (' ' == *s)
    {
#ifdef DEBUG_BUILD_FROM_SMARTS
      cerr << "Found space, ending parsing\n";
#endif
      if (0 != paren_level)
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Un-closed parenthesis");
        return 0;
      }

      return 1;
    }
    else if ('(' == *s)
    {
      if (! check_compatiability_table(previous_token_was, PREVIOUS_TOKEN_WAS_OPEN_PAREN))
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, 
                              "Branch specifier can only follow an atom or ring specifier");
        return 0;
      }

      assert(!previous_bond);

//    Push the various stacks.

      paren_level++;

      atom_stack.add(previous_atom);
      chirality_stack.add(previous_atom_chiral_count);

      characters_processed++;
    }
    else if (characters_processed && ')' == *s)
    {
      if (! check_compatiability_table(previous_token_was, PREVIOUS_TOKEN_WAS_CLOSE_PAREN))
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Close paren can only follow an atom");
        return 0;
      }
      if (atom_stack.empty())
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Parenthesis mismatch");
        return 0;
      }

      if (0 == paren_level)
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Too many close paren's");
        return 0;
      }

//    Jun 98. Handle something like '[CD3](:0)'   note it is a digit '0' rather than 'O' (oxygen)

      if (previous_bond)
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Closing paren after bond??");
        return 0;
      }

      paren_level--;

      previous_atom = atom_stack.pop();
      previous_token_was = PREVIOUS_TOKEN_WAS_ATOM;
      previous_atom_chiral_count = chirality_stack.pop();
      characters_processed++;
    }
    else if ('.' == *s)
    {
      smiles_error_message(smarts, characters_to_process, characters_processed, 
                            "disconnect '.' not allowed in smarts");
      return 0;
    }
    else if (characters_processed && ('/' == *s))
    {
      if (! check_compatiability_table(previous_token_was, PREVIOUS_TOKEN_WAS_BOND))
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Incorrectly placed / bond");
        return 0;
      }
       
      characters_processed++;

      previous_bond.reset(new Substructure_Bond());
      previous_bond->make_single_or_aromatic();
    }
    else if (characters_processed && ('\\' == *s))
    {
      if (! check_compatiability_table(previous_token_was, PREVIOUS_TOKEN_WAS_BOND))
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Incorrectly placed \\ bond");
        return 0;
      }
       
      characters_processed++;

      previous_bond.reset(new Substructure_Bond());
      previous_bond->make_single_or_aromatic();
    }
    else
    {
      if (! check_compatiability_table(previous_token_was, PREVIOUS_TOKEN_WAS_ATOM))
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "incorrectly placed atom");
        return 0;
      }

      Substructure_Atom * a;
      if (0 == atoms_this_fragment)
        a = this;
      else
        a = new Substructure_Atom;

      int tmp;
      if ('[' == *s)
        tmp = a->construct_from_smarts_token(s, characters_to_process - characters_processed);
      else if ('*' == *s)
        tmp = 1;
      else
        tmp = a->construct_from_smiles_token(s, characters_to_process - characters_processed);

#ifdef DEBUG_BUILD_FROM_SMARTS
      cerr << "After processing '" << s << "' tmp = " << tmp << '\n';
#endif

      if (0 == tmp)
      {
        smiles_error_message(smarts, characters_to_process, characters_processed, "Cannot parse smarts");
        if (atoms_this_fragment > 0)      // only delete A if we created it
          delete a;
        return 0;
      }

      completed[atoms_in_previous_disconnected_sections + atoms_this_fragment] = a;
      pst.set_last_query_atom_created(atoms_in_previous_disconnected_sections + atoms_this_fragment);

      if (a->initial_atom_number() < 0)    // only change it if it is unset
        a->set_initial_atom_number(atoms_in_previous_disconnected_sections + atoms_this_fragment);

      a->set_unique_id (atoms_in_previous_disconnected_sections + atoms_this_fragment);

      atoms_this_fragment++;

      if (nullptr != previous_atom)
      {
        Substructure_Bond * b;

        if (previous_bond)
          b = previous_bond.release();
        else
        {
          b = new Substructure_Bond();
          b->make_single_or_aromatic();
        }

#ifdef DEBUG_BUILD_FROM_SMARTS
        cerr << "Making bond between atom " << previous_atom->unique_id() << " and " << atoms_this_fragment - 1 << '\n';
#endif

        b->set_atom(previous_atom);
        a->_add_bond(b);
      }

      previous_atom = a;

//    No chirality in smarts

      characters_processed += tmp;
    }
//  cerr << "Character '" << *s << " interpreted as " << previous_token_was << " nchars = " << characters_processed << '\n';
  }

  if (0 != paren_level)
  {
    smiles_error_message(smarts, characters_to_process, characters_processed, "mismatched parentheses");
    return 0;
  }

  if (PREVIOUS_TOKEN_WAS_BOND == previous_token_was)
  {
    smiles_error_message(smarts, characters_to_process, characters_processed, "a smarts cannot end with a bond");
    return 0;
  }

#ifdef DEBUG_BUILD_FROM_SMARTS
  cerr << "Substructure_Atom::_parse_smarts_specifier: returning 1\n";
#endif

  return 1;
}

/*
  A single Substructure_Atom object can parse a complete smarts, because
  connected atoms become children
*/

int
Substructure_Atom::parse_smarts_specifier(const const_IWSubstring & smarts,
                                          Parse_Smarts_Tmp & pst,
                                          int atoms_in_previous_disconnected_sections)
{
  const_IWSubstring mysmarts(smarts);
  mysmarts.strip_leading_blanks();

  Smiles_Ring_Status ring_status;

  int rc = _parse_smarts_specifier(mysmarts, pst, atoms_in_previous_disconnected_sections, ring_status);

  _attributes_specified =  Substructure_Atom_Specifier::count_attributes_specified();

  count_attributes_specified();  // Recursively initialise all atoms.

  ok_recursive();

  if (0 == rc)
    return 0;

  if (! ring_status.complete())
  {
    ring_status.report_hanging_ring_closures(cerr);
    return 0;
  }

  return rc;
}

int
Substructure_Atom::parse_smarts_specifier(const const_IWSubstring & smarts)
{
  Parse_Smarts_Tmp pst;

  pst.set_natoms(smarts.length());

  int notused = 0;

  return parse_smarts_specifier(smarts, pst, notused);
}

void
reset_parse_smarts_file_scope_variables() {
  if (nullptr != compatability_table) {
    delete [] compatability_table;
    compatability_table = nullptr;
  }
}
