load("@nanobind_bazel//:build_defs.bzl", "nanobind_extension")
load("@rules_python//python:defs.bzl", "py_test")

package(default_visibility = ["//visibility:public"])

nanobind_extension(
    name = "lillymol_nb",
    srcs = [
        "lillymol_nb.cc",
        "lillymol_nb_atom_bond.cc",
        "lillymol_nb_chirality.cc",
        "lillymol_nb_common.cc",
        "lillymol_nb_descriptors.cc",
        "lillymol_nb_fingerprint.cc",
        "lillymol_nb_internal.h",
        "lillymol_nb_io.cc",
        "lillymol_nb_molecule.cc",
        "lillymol_nb_set_of_atoms.cc",
        "lillymol_nb_standardise.cc",
        "lillymol_nb_substructure.cc",
        "lillymol_nb_tools.cc",
        "lillymol_nb_tsubstructure.cc",
    ],
    deps = [
        "//Molecule_Lib:atom_pair_fingerprint",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:linear_fingerprint",
        "//Molecule_Lib:moleculeio",
        "//pybind:tsubstructure",
        "//Molecule_Tools:alogp_lib",
        "//Molecule_Tools:dicer_api",
        "//Molecule_Tools:ec_fingerprint_lib",
        "//Molecule_Tools:jwcats_lib",
        "//Molecule_Tools:nvrtspsa",
        "//Molecule_Tools:qed_lib",
        "//Molecule_Tools:unique_molecules_api",
        "//Molecule_Tools:xlogp_lib",
        "//Utilities/GFP_Tools:truncated_distance_matrix",
    ],
)

py_test(
    name = "lillymol_nb_test",
    srcs = ["lillymol_nb_test.py"],
    data = [":lillymol_nb"],
    deps = ["//Utilities/GFP_Tools:nearneighbours_py_proto"],
)
