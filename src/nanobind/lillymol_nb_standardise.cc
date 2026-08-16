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

  nb::class_<Element_Transformations>(m, "ElementTransformations")
      .def(nb::init<>())
      .def("active", [](const Element_Transformations& etrans) {
        return static_cast<bool>(etrans.active());
      }, "True if any transformations have been added")
      .def("add",
           [](Element_Transformations& etrans, const std::string& directive) {
             const IWString tmp(directive);
             return static_cast<bool>(etrans.Add(tmp));
           },
           nb::arg("directive"), "Add a transformation directive like 'I=Cl'")
      .def("process",
           [](Element_Transformations& etrans, Molecule& mol) { return etrans.process(mol); },
           nb::arg("mol"), "Apply transformations to molecule");
}

}  // namespace lillymol_nb
