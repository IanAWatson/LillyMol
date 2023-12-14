cc_library (
    name = "lillymol_julia",
    srcs = [
        "lillymol.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@libcxxwrap//:libcxxwrap",
        "@julia//:julia",
    ],
)

cc_library(
    name = "bag_julia",
    srcs = [
        "bag.cc",
    ],
    deps = [
        "@libcxxwrap//:libcxxwrap",
        "@julia//:julia",
    ],
)

cc_shared_library(
    name = "lillymol_shared",
    deps = [
        "lillymol_julia",
    ],
    dynamic_deps = [
    ],
    shared_lib_name = "lillymol_julia.so",
)

cc_shared_library(
    name = "box_shared",
    deps = [
        ":bag_julia",
    ],
    dynamic_deps = [
    ],
    shared_lib_name = "bag_julia.so",
)
