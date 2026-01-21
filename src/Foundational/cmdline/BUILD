load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")
cc_library(
    name = "iwcmdline",
    srcs = [
        "cmdline.cc",
        "cmdline_stdstring.cc",
    ],
    hdrs = [
        "cmdline.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//Foundational/iwmisc:iwconfig",
        "//Foundational/iwstring",
    ],
)
