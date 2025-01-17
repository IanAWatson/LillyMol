load("@com_google_protobuf//:protobuf.bzl", "py_proto_library")
load("@rules_proto//proto:defs.bzl", "proto_library")
# load("@rules_proto_grpc_go//:defs.bzl", "go_proto_compile", "go_proto_library")
load("//build_deps:install.bzl", "local_install")

local_install(
    name = "install",
    srcs = [
        ":xgboost_model_evaluate",
    ],
)

proto_library(
    name = "xgboost_model_proto",
    srcs = [
        "xgboost_model.proto",
    ],
    tags = [
        "xgboost",
    ],
)

cc_proto_library(
    name = "xgboost_model_cc_proto",
    tags = [
        "xgboost",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        ":xgboost_model_proto",
    ],
)

py_proto_library(
    name = "xgboost_model_py_proto",
    srcs = [
        "xgboost_model.proto",
    ],
    tags = [
        "xgboost",
    ],
)

#go_proto_compile(
#    name = "xgboost_model_go_proto",
#    protos = [
#        ":xgboost_model_proto",
#    ],
#)

#go_proto_library(
#    name = "xgboost_model_go_proto",
#    importpath = "xgboost_lib",
#    protos = [
#        ":xgboost_model_proto",
#    ],
#    visibility = ["//visibility:public"],
#)

proto_library(
    name = "random_forest_model_proto",
    srcs = [
        "random_forest_model.proto",
    ],
)

cc_proto_library(
    name = "random_forest_model_cc_proto",
    deps = [
        "random_forest_model_proto",
    ],
)

py_proto_library(
    name = "random_forest_model_py_proto",
    srcs = [
        "random_forest_model.proto",
    ],
)

cc_binary(
    name = "xgboost_model_evaluate",
    srcs = [
        "xgboost_model_evaluate.cc",
    ],
    tags = [
        "xgboost",
    ],
    deps = [
        ":xgboost_model_cc_proto",
        "//Foundational/cmdline_v2",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "@third_party//:dlmc",
        "@third_party//:xgboost",
    ],
)
