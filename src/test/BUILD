load("@com_google_protobuf//bazel:py_proto_library.bzl", "py_proto_library")
load("@rules_proto//proto:defs.bzl", "proto_library")
load("@rules_python//python:defs.bzl", "py_test")

package(default_visibility = ["//visibility:public"])

proto_library(
    name = "lillymol_tests_proto",
    srcs = ["lillymol_tests.proto"],
)

py_proto_library(
    name = "lillymol_tests_py_proto",
    deps = [":lillymol_tests_proto"],
)

# The executable-test runner uses Ruby outside Bazel. Generate its bindings
# with the same compiler Bazel uses for the Python and C++ protos.
genrule(
    name = "lillymol_tests_ruby_proto",
    srcs = ["lillymol_tests.proto"],
    outs = ["lillymol_tests_pb.rb"],
    cmd = "$(execpath @com_google_protobuf//:protoc) " +
          "--proto_path=test --ruby_out=$(@D) $(location lillymol_tests.proto)",
    tools = ["@com_google_protobuf//:protoc"],
)

filegroup(
    name = "generated_protos",
    srcs = [
        ":lillymol_tests_py_proto",
        ":lillymol_tests_ruby_proto",
    ],
)

py_test(
    name = "lillymol_tests_proto_test",
    size = "small",
    srcs = ["lillymol_tests_proto_test.py"],
    deps = [":lillymol_tests_py_proto"],
)
