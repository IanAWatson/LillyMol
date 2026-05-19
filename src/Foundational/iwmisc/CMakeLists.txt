add_library(iwconfig INTERFACE)

target_sources(iwconfig
  INTERFACE
    iwconfig.h
)

target_include_directories(iwconfig
  INTERFACE
    ${LILLYMOL_SRC_DIR}
)

add_library(iwmisc_internal INTERFACE)

target_sources(iwmisc_internal
  INTERFACE
    new_array_.h
)

target_include_directories(iwmisc_internal
  INTERFACE
    ${LILLYMOL_SRC_DIR}
)

if(TARGET tomlplusplus::tomlplusplus)
  add_library(ltoml INTERFACE)

  target_sources(ltoml
    INTERFACE
      toml_support.h
  )

  target_include_directories(ltoml
    INTERFACE
      ${LILLYMOL_SRC_DIR}
  )

  target_link_libraries(ltoml
    INTERFACE
      absl::status
      tomlplusplus::tomlplusplus
  )
else()
  add_library(ltoml INTERFACE)
endif()

add_library(compile_time
  compile_time.cc
)

target_include_directories(compile_time
  PUBLIC
    ${LILLYMOL_SRC_DIR}
)

add_library(activity_from_file INTERFACE)

target_sources(activity_from_file
  INTERFACE
    activity_data_from_file.h
)

target_link_libraries(activity_from_file
  INTERFACE
    iwcmdline
    iwstring_data_source
    iwstring
)

add_library(iwmisc
  _iw_write_array_int.cc
  _iwmaxid_double_int.cc
  _iwmaxid_float_iwstring.cc
  _iwmaxid_int_int.cc
  _iwminid_double_int.cc
  _iwminmax_float.cc
  _iwqsort.cc
  _numeric_data_from_file_float.cc
  _numeric_data_from_file_int.cc
  _set_or_unset_float.cc
  append_first_token_of_name.cc
  dash_d.cc
  dash_f.cc
  dash_s.cc
  dash_x.cc
  endian.cc
  fetch_via_regexp.cc
  fname.cc
  fraction_as_string.cc
  hton_stuff.cc
  int_comparator.cc
  iwdigits.cc
  iwdirname.cc
  logical_expression.cc
  msi_object.cc
  new_double.cc
  new_float.cc
  new_int.cc
  new_unsigned_int.cc
  proto_support.cc
  report_progress.cc
  write_first_token_of_string.cc
  write_space_suppressed_string.cc
)

target_include_directories(iwmisc
  PUBLIC
    ${LILLYMOL_SRC_DIR}
)

target_link_libraries(iwmisc
  PUBLIC
    iwmisc_internal
    iwconfig
    ltoml
    iwcmdline
    iwstring_data_source
    iwqsort
    protobuf::libprotobuf
)

add_library(iwre2
  iwre2.cc
)

target_include_directories(iwre2
  PUBLIC
    ${LILLYMOL_SRC_DIR}
)

target_link_libraries(iwre2
  PUBLIC
    iwstring
    re2::re2
)

add_library(sparse_fp_creator
  sparse_fp_creator.cc
)

target_include_directories(sparse_fp_creator
  PUBLIC
    ${LILLYMOL_SRC_DIR}
)

target_link_libraries(sparse_fp_creator
  PUBLIC
    iwmisc
    iwstring
)

protobuf_generate_cpp(
  NORMALISATION_PROTO_SRCS
  NORMALISATION_PROTO_HDRS
  normalisation.proto
)

add_library(normalisation_cc_proto
  ${NORMALISATION_PROTO_SRCS}
  ${NORMALISATION_PROTO_HDRS}
)

target_include_directories(normalisation_cc_proto
  PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}
)

target_link_libraries(normalisation_cc_proto
  PUBLIC
    protobuf::libprotobuf
)

add_library(normalisation
  normalisation.cc
)

target_include_directories(normalisation
  PUBLIC
    ${LILLYMOL_SRC_DIR}
)

target_link_libraries(normalisation
  PUBLIC
    normalisation_cc_proto
    iwmisc
    iwstring
)

add_library(iw_tabular_data
  iw_tabular_data.cc
)

target_sources(iw_tabular_data
  PUBLIC
    iw_tabular_data.h
)

target_include_directories(iw_tabular_data
  PUBLIC
    ${LILLYMOL_SRC_DIR}
)

target_link_libraries(iw_tabular_data
  PUBLIC
    iwstring_data_source
    iwbits
)

target_link_libraries(normalisation
  PUBLIC
    normalisation_cc_proto
    iwmisc
    iwstring
)
