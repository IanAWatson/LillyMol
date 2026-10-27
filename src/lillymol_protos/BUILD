load("@rules_python//python:defs.bzl", "py_library", "py_test")

package(default_visibility = ["//visibility:public"])

# Put the generated files at their Python package paths for installation.
# Each input was generated with import_prefix="lillymol_protos", so message
# __module__ names also refer to the installed package (including for pickle).
genrule(
    name = "python_protos",
    srcs = [
        "package_init.py",
        "//xgboost:xgboost_model_script_py_proto",
        "//xgboost:random_forest_model_script_py_proto",
        "//Utilities/General:class_label_translation_script_py_proto",
        "//Utilities/General:linear_scaling_script_py_proto",
    ],
    outs = [
        "__init__.py",
        "xgboost_model_pb2.py",
        "random_forest_model_pb2.py",
        "class_label_translation_pb2.py",
        "linear_scaling_pb2.py",
    ],
    cmd = "cp $(location package_init.py) $(@D)/__init__.py && " +
          "cp $(location //xgboost:xgboost_model_script_py_proto) $(@D)/xgboost_model_pb2.py && " +
          "cp $(location //xgboost:random_forest_model_script_py_proto) $(@D)/random_forest_model_pb2.py && " +
          "cp $(location //Utilities/General:class_label_translation_script_py_proto) $(@D)/class_label_translation_pb2.py && " +
          "cp $(location //Utilities/General:linear_scaling_script_py_proto) $(@D)/linear_scaling_pb2.py",
)

py_library(
    name = "python",
    srcs = [":python_protos", "check_runtime.py"],
    deps = [
        "//xgboost:xgboost_model_script_py_proto",
        "//xgboost:random_forest_model_script_py_proto",
        "//Utilities/General:class_label_translation_script_py_proto",
        "//Utilities/General:linear_scaling_script_py_proto",
    ],
)

filegroup(
    name = "generated_protos",
    srcs = [":python_protos", "//xgboost:model_ruby_protos"],
)

py_test(
    name = "metadata_test",
    size = "small",
    srcs = ["metadata_test.py"],
    data = ["testdata/legacy_metadata.json"],
    deps = [":python"],
)
