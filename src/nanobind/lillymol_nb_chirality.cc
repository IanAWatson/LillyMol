#include "nanobind/lillymol_nb_internal.h"

namespace lillymol_nb {

void
BindChirality(nb::module_& m) {
  nb::enum_<ChiralType>(m, "ChiralType")
      .value("CHI_UNSPECIFIED", ChiralType::kChiUnspecified)
      .value("CHI_TETRAHEDRAL_CW", ChiralType::kChiTetrahedralCw)
      .value("CHI_TETRAHEDRAL_CCW", ChiralType::kChiTetrahedralCcw)
      .value("CHI_OTHER", ChiralType::kChiOther);

  nb::enum_<CahnIngoldPrelog>(m, "CIP")
      .value("R", CahnIngoldPrelog::R)
      .value("S", CahnIngoldPrelog::S)
      .value("Neither", CahnIngoldPrelog::kNeither)
      .value("Unspecified", CahnIngoldPrelog::kUnspecified);

  nb::class_<Chiral_Centre>(m, "Chiral_Centre")
      .def(nb::init<atom_number_t>(), nb::arg("atom"))
      .def("atom", &Chiral_Centre::a)
      .def("top_front", &Chiral_Centre::top_front)
      .def("top_back", &Chiral_Centre::top_back)
      .def("left_down", &Chiral_Centre::left_down)
      .def("right_down", &Chiral_Centre::right_down)
      .def("invert", &Chiral_Centre::invert)
      .def("involves",
           [](const Chiral_Centre& centre, atom_number_t atom) {
             return static_cast<bool>(centre.involves(atom));
           },
           nb::arg("atom"))
      .def("implicit_hydrogens", &Chiral_Centre::implicit_hydrogen_count)
      .def("lone_pairs", &Chiral_Centre::lone_pair_count)
      .def("implicit_hydrogen_is_now_atom_number",
           &Chiral_Centre::implicit_hydrogen_is_now_atom_number)
      .def("lone_pair_is_now_atom_number", &Chiral_Centre::lone_pair_is_now_atom_number)
      .def("atom_is_now_implicit_hydrogen", &Chiral_Centre::atom_is_now_implicit_hydrogen)
      .def("atom_is_now_lone_pair", &Chiral_Centre::atom_is_now_lone_pair)
      .def("atom_numbers_are_swapped", &Chiral_Centre::atom_numbers_are_swapped)
      .def("__repr__", &ChiralCentreRepr);

  m.def("tetrahedral_chirality", &TetrahedralChirality,
        nb::arg("mol"), nb::arg("atom"), nb::arg("check_is_chiral") = false);
  m.def("is_actually_chiral",
        [](Molecule& mol, atom_number_t atom) {
          if (atom < 0 || atom >= mol.natoms()) {
            throw std::invalid_argument("atom number outside molecule");
          }
          return static_cast<bool>(::is_actually_chiral(mol, atom));
        },
        nb::arg("mol"), nb::arg("atom"));
  m.def("is_chiral_implicit_hydrogen",
        [](int connection) { return IsChiralImplicitHydrogen(connection); },
        nb::arg("connection"));
  m.def("is_chiral_lone_pair",
        [](int connection) { return IsChiralLonePair(connection); },
        nb::arg("connection"));
}

}  // namespace lillymol_nb
