protobuf_generate_cpp(EVIDENCE_PROTO_SRCS EVIDENCE_PROTO_HDRS evidence.proto)
add_library(evidence_cc_proto ${EVIDENCE_PROTO_SRCS} ${EVIDENCE_PROTO_HDRS})
target_include_directories(evidence_cc_proto PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(evidence_cc_proto PUBLIC protobuf::libprotobuf)

protobuf_generate_cpp(GFP_MODEL_PROTO_SRCS GFP_MODEL_PROTO_HDRS gfp_model.proto)
add_library(gfp_model_cc_proto ${GFP_MODEL_PROTO_SRCS} ${GFP_MODEL_PROTO_HDRS})
target_include_directories(gfp_model_cc_proto PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(gfp_model_cc_proto PUBLIC protobuf::libprotobuf)

protobuf_generate_cpp(GFP_TO_SVM_LITE_PROTO_SRCS GFP_TO_SVM_LITE_PROTO_HDRS gfp_to_svm_lite.proto)
add_library(gfp_to_svm_lite_cc_proto ${GFP_TO_SVM_LITE_PROTO_SRCS} ${GFP_TO_SVM_LITE_PROTO_HDRS})
target_include_directories(gfp_to_svm_lite_cc_proto PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(gfp_to_svm_lite_cc_proto PUBLIC protobuf::libprotobuf)

protobuf_generate_cpp(NEARNEIGHBOURS_PROTO_SRCS NEARNEIGHBOURS_PROTO_HDRS nearneighbours.proto)
add_library(nearneighbours_cc_proto ${NEARNEIGHBOURS_PROTO_SRCS} ${NEARNEIGHBOURS_PROTO_HDRS})
target_include_directories(nearneighbours_cc_proto PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(nearneighbours_cc_proto PUBLIC protobuf::libprotobuf)

protobuf_generate_cpp(NN_REQUEST_PROTO_SRCS NN_REQUEST_PROTO_HDRS nn_request.proto)
add_library(nn_request_cc_proto ${NN_REQUEST_PROTO_SRCS} ${NN_REQUEST_PROTO_HDRS})
target_include_directories(nn_request_cc_proto PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(nn_request_cc_proto PUBLIC nearneighbours_cc_proto protobuf::libprotobuf)

add_library(gfp
  _fixed_size_counted_fingerprint_uchar.cc
  _fixed_size_counted_fingerprint_uint.cc
  bic_table.cc
  bit_subset.cc
  count_tdts.cc
  dyfp.cc
  fixed_size_counted_fingerprint_uchar.cc
  gfp.cc
  multi_conformer.cc
  nndata.cc
  sparse_collection.cc
  sparsefp.cc
  tversky.cc
  various_distance_metrics.cc
)

target_include_directories(gfp PUBLIC ${LILLYMOL_SRC_DIR})

target_link_libraries(gfp
  PUBLIC
    gfp_to_svm_lite_cc_proto
    nearneighbours_cc_proto
    accumulator
    iwcmdline
    iwstring_data_source
    iw_tdt
    iwbits
    iwmisc
    sparse_fp_creator
    iwqsort
    re2::re2
)

add_library(smiles_id_dist
  smiles_id_dist.cc
)

target_include_directories(smiles_id_dist PUBLIC ${LILLYMOL_SRC_DIR})
target_link_libraries(smiles_id_dist PUBLIC gfp)

add_library(gfp_standard
  gfp_standard.cc
)

target_include_directories(gfp_standard PUBLIC ${LILLYMOL_SRC_DIR})
target_link_libraries(gfp_standard PUBLIC gfp smiles_id_dist iwmisc)

add_library(gfp_bit_subset
  gfp_bit_subset.cc
)

target_include_directories(gfp_bit_subset PUBLIC ${LILLYMOL_SRC_DIR})
target_link_libraries(gfp_bit_subset PUBLIC gfp)

add_library(neighbour_list
  _neighbour_list_fp_and_smiles.cc
  neighbour_list.cc
)

target_include_directories(neighbour_list PUBLIC ${LILLYMOL_SRC_DIR})
target_link_libraries(neighbour_list PUBLIC smiles_id_dist)

add_library(spread spread.cc)
target_include_directories(spread PUBLIC ${LILLYMOL_SRC_DIR})
target_link_libraries(spread PUBLIC gfp smiles_id_dist)

add_library(build_pool build_pool.cc)
target_include_directories(build_pool PUBLIC ${LILLYMOL_SRC_DIR})
target_link_libraries(build_pool PUBLIC gfp smiles_id_dist)

add_library(distance_scaling distance_scaling.cc)
target_include_directories(distance_scaling PUBLIC ${LILLYMOL_SRC_DIR})
target_link_libraries(distance_scaling PUBLIC iwmisc)

add_library(bsquared bsquared.cc)
target_include_directories(bsquared PUBLIC ${LILLYMOL_SRC_DIR})
target_link_libraries(bsquared PUBLIC accumulator)

add_library(metric
  Enrichment.cc
  Metric.cc
)
target_include_directories(metric PUBLIC ${LILLYMOL_SRC_DIR})
target_link_libraries(metric PUBLIC accumulator)

add_library(bit_and_weight bit_and_weight.cc)
target_include_directories(bit_and_weight PUBLIC ${LILLYMOL_SRC_DIR})
target_link_libraries(bit_and_weight PUBLIC gfp iwstring)

add_library(nbc nbc.cc)
target_include_directories(nbc PUBLIC ${LILLYMOL_SRC_DIR})
target_link_libraries(nbc PUBLIC bit_and_weight gfp iwstring_data_source iwstring)

add_library(bit_subset_v2 INTERFACE)
target_sources(bit_subset_v2 INTERFACE bit_subset_v2.h)
target_include_directories(bit_subset_v2 INTERFACE ${LILLYMOL_SRC_DIR})
target_link_libraries(bit_subset_v2 INTERFACE gfp iwstring_data_source)

add_library(spread_weights spread_weights.cc)
target_include_directories(spread_weights PUBLIC ${LILLYMOL_SRC_DIR})
target_link_libraries(spread_weights PUBLIC iwcmdline iwstring_data_source iwstring)

if(LILLYMOL_BUILD_TESTS)
  add_executable(sparsefp_test sparsefp_test.cc)
  target_link_libraries(sparsefp_test PRIVATE gfp GTest::gtest_main)
  add_test(NAME sparsefp_test COMMAND sparsefp_test)
endif()


function(lillymol_gfp_executable name)
  add_executable(${name} ${ARGN})
  target_include_directories(${name} PRIVATE ${LILLYMOL_SRC_DIR})
endfunction()

lillymol_gfp_executable(gfp_distance_matrix gfp_distance_matrix.cc)
target_link_libraries(gfp_distance_matrix PRIVATE
  gfp gfp_standard iwcmdline iwstring_data_source iwmisc
)

lillymol_gfp_executable(gfp_nearneighbours gfp_nearneighbours.cc)
target_link_libraries(gfp_nearneighbours PRIVATE
  gfp gfp_standard neighbour_list iwcmdline iwstring_data_source iwmisc
)

lillymol_gfp_executable(gfp_spread gfp_spread.cc)
target_link_libraries(gfp_spread PRIVATE
  gfp gfp_standard spread build_pool iwcmdline iwstring_data_source iwmisc
)

lillymol_gfp_executable(gfp_leader gfp_leader.cc)
target_link_libraries(gfp_leader PRIVATE
  gfp gfp_standard iwcmdline iwstring_data_source iwmisc
)

lillymol_gfp_executable(gfp_to_svm_lite gfp_to_svm_lite.cc)
target_link_libraries(gfp_to_svm_lite PRIVATE
  gfp gfp_standard gfp_to_svm_lite_cc_proto iwcmdline iwstring_data_source iwmisc
)
