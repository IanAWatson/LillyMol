add_library(iwbits
  bits_in_common.c
  bits_in_common_avx.cc
  du_bin2ascii.cc
  fastprint.cc
  fixed_bit_vector.cc
  fixed_size_counted_fingerprint.cc
  hexbits.cc
  iwbits.cc
  swap_byte_order.cc
)

target_include_directories(iwbits
  PUBLIC
    ${LILLYMOL_SRC_DIR}
)

target_compile_options(iwbits
  PRIVATE
    -mavx512f
    -mavx512vpopcntdq
)

target_link_libraries(iwbits
  PUBLIC
    iwstring_data_source
    iwmisc
    iwconfig
    iwstring
)

if(LILLYMOL_BUILD_BENCHMARKS)
  add_executable(fixed_bit_vector_benchmark
    fixed_bit_vector_benchmark.cc
  )

  target_link_libraries(fixed_bit_vector_benchmark
    PRIVATE
      iwbits
      benchmark::benchmark
  )
endif()

if(LILLYMOL_BUILD_TESTS)
  add_executable(fixed_bit_vector_test
    fixed_bit_vector_test.cc
  )

  target_link_libraries(fixed_bit_vector_test
    PRIVATE
      iwbits
      GTest::gtest
      GTest::gtest_main
  )

  add_test(NAME fixed_bit_vector_test COMMAND fixed_bit_vector_test)
endif()
