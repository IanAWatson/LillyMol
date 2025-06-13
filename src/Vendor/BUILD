
cc_binary(
    name = "clogp2descriptors_biobyte",
    srcs = [
        "clogp2descriptors_biobyte.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "iwclogp",
    srcs = [
        "biobyte_stuff.c",
        "iwclogp.cc",
    ],
    tags = [
        "vendor",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@biobyte",
    ],
)

cc_binary(
    name = "marvin_pka",
    srcs = [
        "marvin_pka.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "rcorina",
    srcs = [
        "rcorina.cc",
    ],
    tags = [
        "vendor",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/string_data_source",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@corina",
    ],
)
