#include "nanobind/lillymol_nb_internal.h"

#include "Molecule_Lib/standardise.h"

namespace lillymol_nb {

void
BindStandardise(nb::module_& m) {
  nb::class_<Chemical_Standardisation>(m, "Standardise")
      .def(nb::init<>())
      .def("activate_all", &Chemical_Standardisation::activate_all,
           "Activate all transformations")
      .def("set_verbose", &Chemical_Standardisation::set_verbose,
           "Set verbosity")
      .def("process", &Chemical_Standardisation::process,
           "Apply active transformations to molecule");
}

}  // namespace lillymol_nb
