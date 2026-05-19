add_library(iwstring_data_source
  iwmmap.cc
  iwstring_data_source.cc
  string_data_source.cc
)

target_include_directories(iwstring_data_source
  PUBLIC
    ${LILLYMOL_SRC_DIR}
)

target_link_libraries(iwstring_data_source
  PUBLIC
    iwaray
    iwconfig
    iwre2
    iwstring
    re2::re2
    ZLIB::ZLIB
)

if(LILLYMOL_ABSL_CRC32C)
  add_library(iwtfdata_record
    tfdatarecord.cc
  )

  target_include_directories(iwtfdata_record
    PUBLIC
      ${LILLYMOL_SRC_DIR}
  )

  target_link_libraries(iwtfdata_record
    PUBLIC
      iwmisc
      iwstring
      ${LILLYMOL_ABSL_CRC32C}
      ${LILLYMOL_PROTOBUF_LIB}
  )
else()
  message(STATUS "Skipping iwtfdata_record")
endif()

add_library(iwrecordio
  iwrecordio.cc
)

target_include_directories(iwrecordio
  PUBLIC
    ${LILLYMOL_SRC_DIR}
)

target_link_libraries(iwrecordio
  PUBLIC
    iwmisc
    iwstring
    protobuf::libprotobuf
)

protobuf_generate_cpp(
  PROTO_FOR_TESTING_SRCS
  PROTO_FOR_TESTING_HDRS
  proto_for_testing.proto
)

add_library(proto_for_testing_cc_proto
  ${PROTO_FOR_TESTING_SRCS}
  ${PROTO_FOR_TESTING_HDRS}
)

target_include_directories(proto_for_testing_cc_proto
  PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}
)

target_link_libraries(proto_for_testing_cc_proto
  PUBLIC
    protobuf::libprotobuf
)

if(LILLYMOL_BUILD_TESTS)
  add_executable(iwstring_data_source_test
    iwstring_data_source_test.cc
  )

  target_link_libraries(iwstring_data_source_test
    PRIVATE
      iwstring_data_source
      GTest::gtest_main
  )

  add_test(NAME iwstring_data_source_test COMMAND iwstring_data_source_test)

  add_executable(iwrecordio_test
    iwrecordio_test.cc
  )

  target_link_libraries(iwrecordio_test
    PRIVATE
      iwrecordio
      proto_for_testing_cc_proto
      GTest::gtest_main
  )

  add_test(NAME iwrecordio_test COMMAND iwrecordio_test)
endif()
