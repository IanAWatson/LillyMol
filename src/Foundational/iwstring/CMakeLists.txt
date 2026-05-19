add_library(iwstring
  IWString_class.cc
  append_number.cc
  hexprint.cc
  iw_stl_hash_map.cc
  iwgetline.cc
  iwstring.cc
  iwstring_and_file_descriptor.cc
  iwstring_strstream.cc
  iwstrncasecmp.cc
  iwwrite.cc
  iwwrite_block.cc
  iwzlib.cc
  remove_suffix.cc
  string_change.cc
  string_relationals.cc
  string_to_float.cc
  tokenise_with_quotes.cc
  unhtml.cc
)

target_include_directories(iwstring
  PUBLIC
    ${LILLYMOL_SRC_DIR}
)

target_link_libraries(iwstring
  PUBLIC
    iwaray
    iwconfig
    ZLIB::ZLIB
    m
)

add_library(absl_hash INTERFACE)

target_sources(absl_hash
  INTERFACE
    absl_hash.h
)

target_link_libraries(absl_hash
  INTERFACE
    iwstring
    absl::hash
)
if(LILLYMOL_BUILD_TESTS)
  add_executable(iwstring_test iwstring_test.cc)
  target_link_libraries(iwstring_test PRIVATE iwstring GTest::gtest_main)
  add_test(NAME iwstring_test COMMAND iwstring_test)

  add_executable(tokenise_with_quotes_test tokenise_with_quotes_test.cc)
  target_link_libraries(tokenise_with_quotes_test PRIVATE iwstring GTest::gtest_main)
  add_test(NAME tokenise_with_quotes_test COMMAND tokenise_with_quotes_test)
endif()
