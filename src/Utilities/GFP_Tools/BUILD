# load("@com_google_protobuf//:protobuf.bzl", "py_proto_library")
# load("@rules_python//python:proto.bzl", "py_proto_library")
load("@com_google_protobuf//bazel:py_proto_library.bzl", "py_proto_library")
# load("@rules_python//python:proto.bzl", "py_proto_library")
# load("@rules_proto_grpc_cpp//:defs.bzl", "cpp_grpc_compile")

load("@rules_proto//proto:defs.bzl", "proto_library")

proto_library(
    name = "evidence_proto",
    srcs = [
        "evidence.proto",
    ],
)

cc_proto_library(
    name = "evidence_cc_proto",
    deps = [
        ":evidence_proto",
    ],
)

proto_library(
    name = "gfp_model_proto",
    srcs = [
        "gfp_model.proto",
    ],
)

cc_proto_library(
    name = "gfp_model_cc_proto",
    deps = [
        ":gfp_model_proto",
    ],
)

proto_library(
    name = "gfp_to_svm_lite_proto",
    srcs = [
        "gfp_to_svm_lite.proto",
    ],
)

cc_proto_library(
    name = "gfp_to_svm_lite_cc_proto",
    deps = [
        ":gfp_to_svm_lite_proto",
    ],
)

py_proto_library(
    name = "gfp_to_svm_lite_proto_py_proto",
    deps = [
        ":gfp_to_svm_lite_proto",
    ],
)

proto_library(
    name = "nearneighbours_proto",
    srcs = [
        "nearneighbours.proto",
    ],
    visibility = [
        "//visibility:public",
    ],
)

cc_proto_library(
    name = "nearneighbours_cc_proto",
    visibility = [
        "//visibility:public",
    ],
    deps = [
        ":nearneighbours_proto",
    ],
)

py_proto_library(
    name = "nearneighbours_py_proto",
    deps = [
        ":nearneighbours_proto",
    ],
)

proto_library(
    name = "nn_request_proto",
    srcs = [
        "nn_request.proto",
    ],
    deps = [
        "//Utilities/GFP_Tools:nearneighbours_proto",
    ],
)

cc_proto_library(
    name = "nn_request_cc_proto",
    deps = [
        ":nn_request_proto",
    ],
)

py_proto_library(
    name = "nn_request_py_proto",
    deps = [
        ":nn_request_proto",
        ":nearneighbours_proto",
    ],
)

#grpc_proto_library(
#  name = "nn_request_grpc",
#  srcs = [
#      "nn_request.proto",
#  ],
#)

cc_binary(
    name = "descriptor_file_to_01_fingerprints",
    srcs = [
        "descriptor_file_to_01_fingerprints.cc",
    ],
    deps = [
        ":gfp",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "evidence",
    srcs = [
        "evidence.cc",
    ],
    deps = [
        ":evidence_cc_proto",
        ":nearneighbours_cc_proto",
        "//Foundational/cmdline_v2",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/data_source:iwtfdata_record",
        "//Foundational/iwstring:absl_hash",
        "@com_google_absl//absl/container:flat_hash_map",
    ],
)

cc_binary(
    name = "gfp_add_descriptors",
    srcs = [
        "gfp_add_descriptors.cc",
    ],
    deps = [
        ":gfp",
    ],
)

# Catboost support is experimental.
# cc_binary(
#     name = "gfp_catboost_score",
#     tags = ["catboost"],
#     srcs = [
#         "gfp_catboost_score.cc",
#         "Catboost.model.bin.cpp",
#     ],
#     deps = [
#         ":gfp",
#         ":gfp_to_svm_lite_cc_proto",
#         ":gfp_model_cc_proto",
#         "//Foundational/cmdline_v2:cmdline_v2",
#         "//Foundational/data_source:iwstring_data_source",
#         "//Utilities/General:class_label_translation_cc_proto",
#         "//Utilities/General:scaler",
#     ],
# )

cc_binary(
    name = "gfp_distance_matrix_iwdm",
    srcs = [
        "gfp_distance_matrix_iwdm.cc",
    ],
    deps = [
        ":gfp",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
        "//Utilities/Distance_Matrix:distance_matrix_base",
    ],
)

cc_binary(
    name = "gfp_distance_filter",
    srcs = [
        "gfp_distance_filter.cc",
    ],
    deps = [
        ":gfp",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
        "//Utilities/Distance_Matrix:distance_matrix_base",
    ],
)

cc_binary(
    name = "gfp_distance_filter_standard",
    srcs = [
        "gfp_distance_filter_standard.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
        "//Utilities/GFP_Tools:gfp_standard",
    ],
)

cc_binary(
    name = "gfp_flatten_counted",
    srcs = [
        "gfp_flatten_counted.cc",
    ],
    deps = [
        ":gfp",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "gfp_histogram",
    srcs = [
        "gfp_histogram.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/histogram:iwhistogram",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
        "//Utilities/GFP_Tools:gfp",
    ],
)

cc_binary(
    name = "gfp_incremental_diversity",
    srcs = [
        "gfp_incremental_diversity.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
        "//Utilities/GFP_Tools:gfp",
    ],
)

cc_binary(
    name = "gfp_iterative_expansion",
    srcs = [
        "gfp_iterative_expansion.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
        "//Utilities/GFP_Tools:gfp",
    ],
)

cc_binary(
    name = "gfp_leader",
    srcs = [
        "gfp_leader.cc",
        "leader.h",
    ],
    deps = [
        ":gfp",
        ":nearneighbours_cc_proto",
        "//Foundational/accumulator",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "gfp_leader_standard",
    srcs = [
        "gfp_leader_standard.cc",
    ],
    copts = [
        #       "-fopenmp=libiomp5",  llvm
        "-fopenmp",
    ],
    linkopts = [
        #       "-L/usr/lib/llvm-11/lib/",  # clang
        "-lgomp",
    ],
    tags = [
        "omp",
    ],
    deps = [
        ":gfp",
        ":gfp_standard",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "gfp_leader_tbb",
    srcs = [
        "gfp_leader_tbb.cc",
        "leader_parallel.cc",
        "leader_parallel.h",
    ],
    linkopts = [
        "-ldl",
    ],
    tags = [
        "tbb",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
        "//Utilities/GFP_Tools:gfp",
        "@onetbb//:tbb",
    ],
)

#cc_binary(
#    name = "gfp_leader_threaded_mpi",
#    srcs = [
#        "gfp_leader_threaded_mpi.cc",
#    ],
#    tags = [
#        "mpi",
#    ],
#    deps = [
#        "//Foundational/accumulator",
#        "//Foundational/cmdline:iwcmdline",
#        "//Foundational/data_source:iwstring_data_source",
#        "//Foundational/iwmisc",
#        "//Utilities/GFP_Tools:gfp",
#        "@mpich",
#    ],
#)

cc_binary(
    name = "gfp_mcs",
    srcs = [
        "gfp_mcs.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/histogram:iwhistogram",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
        "//Utilities/GFP_Tools:gfp",
    ],
)

cc_binary(
    name = "gfp_naive_bayesian",
    srcs = [
        "gfp_naive_bayesian.cc",
        #     "bit_subset_v2.h",
    ],
    deps = [
        ":bit_subset_v2",
        ":gfp",
        ":smiles_id_dist",
        ":spread",
        "//Foundational/accumulator",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "gfp_server",
    srcs = [
        "gfp_server.cc",
    ],
    tags = [
        "gfp_server",
        "zeromq",
    ],
    deps = [
        "nn_request_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Tools:maccskeys_fn5_lib",
        "//Molecule_Tools:mpr",
        "//Utilities/GFP_Tools:gfp_standard",
        "//Utilities/GFP_Tools:nearneighbours_cc_proto",
        "@com_google_protobuf//:protobuf",
        "@third_party//:cppzmq",
        "@third_party//:libzmq",
    ],
)

# cc_binary(
#   name = "gfp_server_grpc",
#   srcs = [
#       "gfp_server_grpc.cc",
#   ],
#   deps = [
#        "nn_request_cc_proto",
#        "//Foundational/cmdline:iwcmdline",
#        "//Foundational/data_source:iwstring_data_source",
#        "//Molecule_Lib:iwmfingerprint",
#        "//Molecule_Lib:iwmolecule",
#        "//Molecule_Lib:moleculeio",
#        "//Molecule_Tools:mpr",
#        "//Molecule_Tools:maccskeys_fn5_lib",
#        "//Utilities/GFP_Tools:gfp_standard",
#        "//Utilities/GFP_Tools:nearneighbours_cc_proto",
#        "@com_github_grpc_grpc//:grpc++",
#        "@com_google_protobuf//:protobuf",
#   ],
#   tags = [
#       "gfp_server",
#       "grpc,
#   ],
# )

cc_binary(
    name = "gfp_sparse_to_fixed",
    srcs = [
        "gfp_sparse_to_fixed.cc",
    ],
    deps = [
        ":gfp",
        ":smiles_id_dist",
        "//Foundational/accumulator",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "gfp_spread",
    srcs = [
        "gfp_spread.cc",
    ],
    deps = [
        ":gfp",
        ":smiles_id_dist",
        ":spread",
        "//Foundational/accumulator",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "gfp_spread_buckets",
    srcs = [
        "gfp_spread_buckets.cc",
    ],
    deps = [
        ":gfp",
        ":smiles_id_dist",
        ":spread",
        "//Foundational/accumulator",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "gfp_spread_omp",
    srcs = [
        "gfp_spread_omp.cc",
    ],
    copts = ["-fopenmp"],
    linkopts = [
        "-lgomp",
    ],
    tags = [
        "omp",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
        "//Utilities/GFP_Tools:spread",
    ],
)

cc_binary(
    name = "gfp_spread_standard",
    srcs = [
        "gfp_spread_standard.cc",
    ],
    copts = [
        #       "-fopenmp=libiomp5", llvm
        "-fopenmp",
    ],
    linkopts = [
        #       "-L/usr/lib/llvm-11/lib/",  # clang
        "-lgomp",
    ],
    tags = [
        "omp",
    ],
    deps = [
        ":gfp",
        ":gfp_standard",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "gfp_standalone",
    srcs = [
        "gfp_standalone.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Tools:maccskeys_fn5_lib",
        "//Molecule_Tools:mpr",
        "//Utilities/GFP_Tools:gfp_standard",
    ],
)

cc_binary(
    name = "gfp_nearneighbours",
    srcs = [
        "gfp_nearneighbours.cc",
    ],
    deps = [
        ":gfp",
        "//Foundational/histogram:iwhistogram",
    ],
)

cc_binary(
    name = "gfp_nearneighbours_single_file",
    srcs = [
        "fp_and_smiles.cc",
        "gfp_nearneighbours_single_file.cc",
    ],
    deps = [
        ":gfp",
        "//Foundational/data_source:iwtfdata_record",
        "//Foundational/histogram:iwhistogram",
    ],
)

cc_binary(
    name = "gfp_nearneighbours_single_file_tbb",
    srcs = [
        "gfp_nearneighbours_single_file_tbb.cc",
    ],
    copts = ["-fopenmp"],
    linkopts = [
        "-ldl",
        "-lgomp",
    ],
    tags = [
        "omp",
        "tbb",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/data_source:iwtfdata_record",
        "//Foundational/histogram:iwhistogram",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sorted_list",
        "//Foundational/iwqsort",
        "//Utilities/GFP_Tools:gfp",
        "@onetbb//:tbb",
    ],
)

cc_binary(
    name = "gfp_lnearneighbours",
    srcs = [
        "fp_and_smiles.cc",
        "gfp_lnearneighbours.cc",
    ],
    deps = [
        ":gfp",
        ":neighbour_list",
    ],
)

cc_binary(
    name = "gfp_lnearneighbours_standard",
    srcs = [
        "fp_and_smiles.cc",
        "gfp_lnearneighbours_standard.cc",
    ],
    copts = [
        #       "-fopenmp=libiomp5",  llvm
        "-fopenmp",
    ],
    linkopts = [
        #       "-L/usr/lib/llvm-11/lib/",  # clang
        "-lgomp",
    ],
    tags = [
        "omp",
    ],
    deps = [
        ":gfp",
        ":gfp_standard",
        ":nearneighbours_cc_proto",
        ":neighbour_list",
        "//Foundational/data_source:iwtfdata_record",
    ],
)

cc_binary(
    name = "gfp_pairwise_distances",
    srcs = [
        "gfp_pairwise_distances.cc",
    ],
    deps = [
        ":gfp",
        "//Foundational/histogram:iwhistogram",
    ],
)

cc_binary(
    name = "gfp_profile_activity_by_bits",
    srcs = [
        "gfp_profile_activity_by_bits.cc",
    ],
    deps = [
        ":gfp",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "gfp_single_linkage",
    srcs = [
        "gfp_single_linkage.cc",
    ],
    deps = [
        ":gfp",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "gfp_svmfp_score",
    srcs = [
        "gfp_svmfp_score.cc",
    ],
    deps = [
        ":gfp_bit_subset",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:normalisation",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Foundational/iwqsort",
        "//Utilities/GFP_Tools:gfp",
    ],
)

cc_binary(
    name = "gfp_svmfp_score_pthread",
    srcs = [
        "gfp_svmfp_score_pthread.cc",
    ],
    deps = [
        ":gfp_bit_subset",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:normalisation",
        "//Utilities/GFP_Tools:gfp",
    ],
)

cc_binary(
    name = "gfp_svmfp_score_tbb",
    srcs = [
        "gfp_svmfp_score_tbb.cc",
    ],
    linkopts = [
        "-ldl",
    ],
    tags = [
        "tbb",
    ],
    deps = [
        "gfp_bit_subset",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:normalisation",
        "//Utilities/GFP_Tools:gfp",
        "@onetbb//:tbb",
    ],
)

cc_binary(
    name = "gfp_svmfp_score_v2",
    srcs = [
        "gfp_svmfp_score_v2.cc",
    ],
    deps = [
        ":gfp_bit_subset",
        ":gfp_model_cc_proto",
        "//Foundational/accumulator",
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:normalisation",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Foundational/iwqsort",
        "//Utilities/General:class_label_translation_cc_proto",
        "//Utilities/General:scaler",
        "//Utilities/GFP_Tools:gfp",
    ],
)

cc_binary(
    name = "gfp_to_descriptors",
    srcs = [
        "gfp_to_descriptors.cc",
    ],
    deps = [
        ":gfp",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "gfp_to_descriptors_multiple",
    srcs = [
        "gfp_to_descriptors_multiple.cc",
    ],
    deps = [
        ":gfp",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "gfp_to_svm_lite.v2",
    srcs = [
        "gfp_to_svm_lite.v2.cc",
    ],
    deps = [
        "gfp_bit_subset",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
        "//Foundational/iwqsort",
        "//Utilities/GFP_Tools:gfp",
    ],
)

cc_binary(
    name = "gfp_to_svm_lite_v3",
    srcs = [
        "gfp_to_svm_lite_v3.cc",
    ],
    deps = [
        ":gfp",
        ":gfp_to_svm_lite_cc_proto",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:activity_from_file",
    ],
)

cc_binary(
    name = "nn_fixed_size_cluster",
    srcs = [
        "nn_fixed_size_cluster.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "nn_merge",
    srcs = [
        "nn_merge.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
        "//Foundational/iwqsort",
        "//Utilities/GFP_Tools:smiles_id_dist",
    ],
)

cc_binary(
    name = "nn_merge_from_smiles",
    srcs = [
        "nn_merge_from_smiles.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
        "//Foundational/iwqsort",
        "//Utilities/GFP_Tools:smiles_id_dist",
    ],
)

cc_binary(
    name = "nn2csv",
    srcs = [
        "nn2csv.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwqsort",
    ],
)

cc_binary(
    name = "nn2proto",
    srcs = [
        "nn2proto.cc",
    ],
    deps = [
        ":gfp",
        ":nearneighbours_cc_proto",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/data_source:iwtfdata_record",
        "//Foundational/iwqsort",
    ],
)

cc_binary(
    name = "nplotnn",
    srcs = [
        "nplotnn.cc",
    ],
    deps = [
        ":distance_scaling",
        ":smiles_id_dist",
        "//Foundational/histogram:iwhistogram",
        "//Foundational/iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "gfp_distance_matrix",
    srcs = [
        "gfp_distance_matrix.cc",
    ],
    deps = [
        ":build_pool",
        ":gfp",
        "//Foundational/histogram:iwhistogram",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "marvin2gfp",
    srcs = [
        "marvin2gfp.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
)

cc_binary(
    name = "nn_leader_and_jp",
    srcs = [
        "nn_leader_and_jp.cc",
    ],
    deps = [
        ":gfp",
        "//Foundational/data_source:iwtfdata_record",
    ],
)

cc_binary(
    name = "nn_spread",
    srcs = [
        "nn_spread.cc",
    ],
    deps = [
        ":nearneighbours_cc_proto",
        ":spread_weights",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwtfdata_record",
    ],
    copts = [
        #       "-fopenmp=libiomp5",  llvm
        "-fopenmp",
    ],
    linkopts = [
        #       "-L/usr/lib/llvm-11/lib/",  # clang
        "-lgomp",
    ],
    tags = [
        "omp",
    ],
)

cc_binary(
    name = "nn_id_to_ndx",
    srcs = [
        "nn_id_to_ndx.cc",
    ],
    deps = [
        ":nearneighbours_cc_proto",
        "//Foundational/accumulator",
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/data_source:iwtfdata_record",
        "//Foundational/iwstring:absl_hash",
    ],
)

cc_binary(
    name = "parallel_nn_search_to_gfp_spread",
    srcs = [
        "parallel_nn_search_to_gfp_spread.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "plate_layout",
    srcs = [
        "plate_layout.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Utilities/GFP_Tools:gfp",
    ],
)

cc_binary(
    name = "random_fingerprint",
    srcs = [
        "random_fingerprint.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwbits",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "test_sparse_bitvector_performance",
    srcs = [
        "test_sparse_bitvector_performance.cc",
    ],
    deps = [
        ":build_pool",
        ":gfp",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "train_test_split_optimise",
    srcs = [
        "train_test_split_optimise.cc",
    ],
    deps = [
        ":nearneighbours_cc_proto",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwtfdata_record",
        "//Foundational/iwmisc",
    ],
)

cc_binary(
    name = "iwstats",
    srcs = [
        "iwstats.cc",
    ],
    deps = [
        ":bsquared",
        ":metric",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc",
        "//Foundational/iwqsort",
    ],
)

cc_binary(
    name = "InformationContent",
    srcs = [
        "InformationContent.cc",
    ],
    deps = [
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwbits",
        "//Foundational/iwmisc",
        "//Utilities/GFP_Tools:gfp",
    ],
)

cc_library(
    name = "nbc",
    srcs = [
        "nbc.cc",
    ],
    hdrs = [
        "nbc.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        ":bit_and_weight",
        ":gfp",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwstring",
    ],
)

cc_library(
    name = "bit_and_weight",
    srcs = [
        "bit_and_weight.cc",
    ],
    hdrs = [
        "bit_and_weight.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        ":gfp",
        "//Foundational/iwstring",
    ],
)

cc_library(
    name = "spread",
    srcs = [
        "spread.cc",
    ],
    hdrs = [
        "spread.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        ":gfp",
        ":smiles_id_dist",
    ],
)

cc_library(
    name = "build_pool",
    srcs = [
        "build_pool.cc",
    ],
    hdrs = [
    ],
    deps = [
        ":gfp",
        ":smiles_id_dist",
    ],
)

cc_library(
    name = "smiles_id_dist",
    srcs = [
        "smiles_id_dist.cc",
    ],
    hdrs = [
        "smiles_id_dist.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        ":gfp",
    ],
)

cc_library(
    name = "gfp",
    srcs = [
        "_fixed_size_counted_fingerprint_uchar.cc",
        "_fixed_size_counted_fingerprint_uint.cc",
        "bic_table.cc",
        "bit_subset.cc",
        "count_tdts.cc",
        "dyfp.cc",
        "fixed_size_counted_fingerprint_uchar.cc",
        "gfp.cc",
        "multi_conformer.cc",
        "nndata.cc",
        "sparse_collection.cc",
        "sparsefp.cc",
        "tversky.cc",
        "various_distance_metrics.cc",
    ],
    hdrs = [
        "bit_subset.h",
        "dyfp.h",
        "fixed_size_counted_fingerprint.h",
        "gfp.h",
        "leader.h",
        "multi_conformer.h",
        "nndata.h",
        "sparse_collection.h",
        "sparsefp.h",
        "tversky.h",
        "various_distance_metrics.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        ":gfp_to_svm_lite_cc_proto",
        ":nearneighbours_cc_proto",
        "//Foundational/accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iw_tdt",
        "//Foundational/iwbits",
        "//Foundational/iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Foundational/iwqsort",
        "@re2",
    ],
)

cc_library(
    name = "gfp_bit_subset",
    srcs = [
        "gfp_bit_subset.cc",
    ],
    hdrs = [
        "gfp_bit_subset.h",
    ],
    deps = [
        ":gfp",
    ],
)

cc_library(
    name = "gfp_standard",
    srcs = [
        "gfp_standard.cc",
    ],
    hdrs = [
        "gfp_standard.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        ":gfp",
        ":smiles_id_dist",
        "//Foundational/iwmisc",
    ],
)

cc_library(
    name = "neighbour_list",
    srcs = [
        "_neighbour_list_fp_and_smiles.cc",
        "neighbour_list.cc",
    ],
    hdrs = [
        "neighbour_list.h",
    ],
    deps = [
        ":smiles_id_dist",
    ],
)

cc_library(
    name = "sparse_collection",
    srcs = [
        "sparse_collection.cc",
    ],
    hdrs = [
        "sparse_collection.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        ":gfp",
        "//Foundational/cmdline:iwcmdline",
    ],
)

cc_library(
    name = "distance_scaling",
    srcs = [
        "distance_scaling.cc",
    ],
    hdrs = [
        "distance_scaling.h",
    ],
    deps = [
        "//Foundational/iwmisc",
    ],
)

cc_library(
    name = "bsquared",
    srcs = [
        "bsquared.cc",
    ],
    hdrs = [
        "bsquared.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        "//Foundational/accumulator",
    ],
)

cc_library(
    name = "bit_subset_v2",
    srcs = [
    ],
    hdrs = [
        "bit_subset_v2.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        ":gfp",
        "//Foundational/data_source:iwstring_data_source",
    ],
)

cc_library(
    name = "metric",
    srcs = [
        "Enrichment.cc",
        "Metric.cc",
    ],
    hdrs = [
        "Enrichment.h",
        "Metric.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        "//Foundational/accumulator",
    ],
)

cc_library( 
    name = "spread_weights",
    srcs = [
        "spread_weights.cc",
    ],
    hdrs = [
        "spread_weights.h",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwstring:iwstring",
    ],
)

cc_binary(
    name = "sparsefp_benchmark",
    srcs = [
        "sparsefp_benchmark.cc",
    ],
    deps = [
        ":gfp",
        "//Foundational/iwmisc:sparse_fp_creator",
        "@benchmark",
    ],
)

#cc_test(
#    name = "bit_subset_test",
#    srcs = [
#        "bit_subset_test.cc",
#    ],
#    timeout="short",
#    deps = [
#        ":gfp",
#        "@googletest//:gtest_main",
#    ],
#)

cc_test(
    name = "sparsefp_test",
    timeout = "short",
    srcs = [
        "sparsefp_test.cc",
    ],
    deps = [
        ":gfp",
        "@googletest//:gtest_main",
    ],
)
