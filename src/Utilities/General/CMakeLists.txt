function(lillymol_general_executable name)
  add_executable(${name} ${ARGN})
  target_include_directories(${name} PRIVATE ${LILLYMOL_SRC_DIR})
endfunction()

protobuf_generate_cpp(CLASS_LABEL_TRANSLATION_SRCS CLASS_LABEL_TRANSLATION_HDRS
  class_label_translation.proto
)

add_library(class_label_translation_cc_proto
  ${CLASS_LABEL_TRANSLATION_SRCS}
  ${CLASS_LABEL_TRANSLATION_HDRS}
)

target_link_libraries(class_label_translation_cc_proto PUBLIC protobuf::libprotobuf)
target_include_directories(class_label_translation_cc_proto PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

protobuf_generate_cpp(FEATURE_SCALING_SRCS FEATURE_SCALING_HDRS
  feature_scaling.proto
)

add_library(feature_scaling_cc_proto
  ${FEATURE_SCALING_SRCS}
  ${FEATURE_SCALING_HDRS}
)

target_link_libraries(feature_scaling_cc_proto PUBLIC protobuf::libprotobuf)
target_include_directories(feature_scaling_cc_proto PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

add_library(scaler INTERFACE)

target_sources(scaler INTERFACE scaler.h)

target_link_libraries(scaler
  INTERFACE
    feature_scaling_cc_proto
)

add_library(linear_scaling
  linear_scaling.cc
)

target_link_libraries(linear_scaling
  PUBLIC
    linear_scaling_cc_proto
    iwmisc
    iwstring
)

add_library(ttitem
  ttitem.cc
)

target_link_libraries(ttitem
  PUBLIC
    iwstring
)

lillymol_general_executable(average
  average.cc iwtokeniser.cc
)
target_link_libraries(average PRIVATE
  accumulator iwstring_data_source iwmisc iwstring
)

lillymol_general_executable(bucketise bucketise.cc)
target_link_libraries(bucketise PRIVATE
  accumulator iwcmdline iwstring_data_source iwmisc
)

lillymol_general_executable(byte_offset_index byte_offset_index.cc)
target_link_libraries(byte_offset_index PRIVATE
  accumulator iwcmdline iwstring_data_source iwmisc
)

lillymol_general_executable(descriptor_file_cat descriptor_file_cat.cc)
target_link_libraries(descriptor_file_cat PRIVATE
  accumulator cmdline_v2 iwstring_data_source iwmisc
)

lillymol_general_executable(descriptor_file_filter descriptor_file_filter.cc)
target_link_libraries(descriptor_file_filter PRIVATE
  accumulator iwcmdline iwstring_data_source iwmisc
)

lillymol_general_executable(descriptor_file_sort descriptor_file_sort.cc)
target_link_libraries(descriptor_file_sort PRIVATE
  accumulator iwcmdline iwstring_data_source iwmisc iwqsort
)

lillymol_general_executable(difference_sort difference_sort.cc)
target_link_libraries(difference_sort PRIVATE
  accumulator cmdline_v2 iwstring_data_source iwmisc iwqsort
)

lillymol_general_executable(distribution distribution.cc)
target_link_libraries(distribution PRIVATE
  accumulator iwcmdline iwstring_data_source iwmisc iwqsort
)

lillymol_general_executable(fetch_sdf fetch_sdf.cc)
target_link_libraries(fetch_sdf PRIVATE
  accumulator iwcmdline iwstring_data_source iwmisc
)

lillymol_general_executable(fetch_sdf_quick fetch_sdf_quick.cc)
target_link_libraries(fetch_sdf_quick PRIVATE
  iwstring_data_source iwmisc iwstring
)

lillymol_general_executable(fetch_smiles
  fetch_smiles.cc iwtokeniser.cc
)
target_link_libraries(fetch_smiles PRIVATE
  iwstring_data_source iwmisc iwstring
)

lillymol_general_executable(fetch_smiles_quick
  fetch_smiles_quick.cc iwtokeniser.cc
)
target_link_libraries(fetch_smiles_quick PRIVATE
  iwstring_data_source iwmisc iwstring
)

lillymol_general_executable(iwsplit iwsplit.cc)
target_link_libraries(iwsplit PRIVATE
  cmdline_v2 iwstring_data_source iwmisc iwstring
)

lillymol_general_executable(iwcut
  iwcut.cc iwtokeniser.cc
)
target_link_libraries(iwcut PRIVATE
  iwstring_data_source iwmisc iwstring
)

lillymol_general_executable(random_column_order random_column_order.cc)
target_link_libraries(random_column_order PRIVATE
  iwcmdline iwstring_data_source iwmisc iwstring
)

lillymol_general_executable(random_records random_records.cc)
target_link_libraries(random_records PRIVATE
  accumulator iwcmdline iwstring_data_source iwmisc iwstring
)

lillymol_general_executable(unique_rows
  unique_rows.cc iwtokeniser.cc
)
target_link_libraries(unique_rows PRIVATE
  accumulator cmdline_v2 iwstring_data_source iwmisc iwqsort iwstring
)

lillymol_general_executable(whatsmissing whatsmissing.cc)
target_link_libraries(whatsmissing PRIVATE
  accumulator cmdline_v2 iwstring_data_source iwmisc iwstring
)

if(LILLYMOL_BUILD_TESTS)
  add_executable(scaler_test scaler_test.cc)

  target_link_libraries(scaler_test PRIVATE
    scaler
    feature_scaling_cc_proto
    GTest::gtest
    GTest::gtest_main
  )

  add_test(NAME scaler_test COMMAND scaler_test)
endif()

lillymol_general_executable(descriptors_to_fingerprint
  descriptors_to_fingerprint.cc
)
target_link_libraries(descriptors_to_fingerprint PRIVATE
  iwstring_data_source
  iwmisc
  sparse_fp_creator
  iwstring
)

lillymol_general_executable(descriptors_to_fingerprint_integer
  descriptors_to_fingerprint_integer.cc
)
target_link_libraries(descriptors_to_fingerprint_integer PRIVATE
  iwcmdline
  iwstring_data_source
  iwmisc
  sparse_fp_creator
  iwstring
  absl_hash
)

lillymol_general_executable(svm_lite_to_gfp
  svm_lite_to_gfp.cc
)
target_link_libraries(svm_lite_to_gfp PRIVATE
  accumulator
  iwcmdline
  iwstring_data_source
  iwmisc
  sparse_fp_creator
)

lillymol_general_executable(iwstats
  iwstats.cc
)
target_link_libraries(iwstats PRIVATE
  accumulator
  iwcmdline
  iwstring_data_source
  iw_tabular_data
  iwmisc
  iwstring
)

lillymol_general_executable(iwstats_summary
  iwstats_summary.cc
)
target_link_libraries(iwstats_summary PRIVATE
  accumulator
  iwcmdline
  iwstring_data_source
  iw_tabular_data
  iwmisc
  iwstring
)

lillymol_general_executable(dfilefilter
  descriptor_file_filter.cc
)
target_link_libraries(dfilefilter PRIVATE
  accumulator
  iwcmdline
  iwstring_data_source
  iw_tabular_data
  iwmisc
  iwstring
)

lillymol_general_executable(correlate
  correlate.cc
)
target_link_libraries(correlate PRIVATE
  accumulator
  iwcmdline
  iwstring_data_source
  iw_tabular_data
  iwmisc
  iwstring
)
