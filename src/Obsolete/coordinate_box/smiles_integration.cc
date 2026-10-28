// Archived SMILES integration excerpts, retained for reference only.
// These functions and options have been removed from Molecule_Lib.
#include "coordinate_box.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/smiles.h"

static int append_coordinate_box_after_each_atom = 0;

void
set_append_coordinate_box_after_each_atom(int s) {
  append_coordinate_box_after_each_atom = s;
}

static void
do_append_coordinate_box(const Atom& a, IWString& smiles) {
  coordinate_box::ConcentricBox box;
  const coordinate_box::LayerPosition layer_position = box.Position(a);
  smiles += "{{B";
  smiles << layer_position;
  smiles += "}}";
}

int
ParseBoxedCoordinates(const_IWSubstring c,  // Local copy.
                      Coordinates& coords) {
  assert(c.starts_with('B'));
  c += 1;  // SKip over 'B'.
  coordinate_box::LayerPosition layer_position;
  const int consumed = coordinate_box::FromString(c, layer_position);
  if (consumed != c.length()) {
    return 0;
  }

  coordinate_box::ConcentricBox box;

  coords = box.CellToCoordinates<float>(layer_position);

  return c.length() + 1;
}

// Former output option branch:
/*
    } else if ("cbox" == tmp) {
      set_append_coordinate_box_after_each_atom(1);
      if (verbose) {
        cerr << "Will add boxed coordinates to each atom in a smiles\n";
      }

*/
