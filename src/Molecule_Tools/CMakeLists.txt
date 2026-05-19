add_library(fileconv_lib
  fileconv_opts.cc
  fix_structures.cc
  known_fragment_data.cc
  remove_duplicate_fragments.cc
)

target_include_directories(fileconv_lib
  PUBLIC
    ${LILLYMOL_SRC_DIR}
)

target_link_libraries(fileconv_lib
  PUBLIC
    moleculeio
    accumulator
    iwcmdline
    iwstring_data_source
    iw_tdt
    iwbits
    iwmisc
    iwstring
)

add_executable(fileconv
  fileconv.cc
  fix_structures.cc
  known_fragment_data.cc
)

target_include_directories(fileconv
  PRIVATE
    ${LILLYMOL_SRC_DIR}
)

target_link_libraries(fileconv
  PRIVATE
    fileconv_lib
    moleculeio
    accumulator
    iwcmdline
    iwstring_data_source
    iw_tdt
    iwbits
    iwmisc
    iwstring
)
