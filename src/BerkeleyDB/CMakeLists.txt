add_library(iwstring_bdb
  iwstring_bdb.cc
)

target_include_directories(iwstring_bdb
  PUBLIC
    ${LILLYMOL_SRC_DIR}
)

target_link_libraries(iwstring_bdb
  PUBLIC
    iwstring
    berkeleydb
)

add_library(unpack_data
  unpack_data.cc
)

target_include_directories(unpack_data
  PUBLIC
    ${LILLYMOL_SRC_DIR}
)

target_link_libraries(unpack_data
  PUBLIC
    iwstring
)

function(lillymol_bdb_executable name)
  add_executable(${name} ${name}.cc)
  target_include_directories(${name} PRIVATE ${LILLYMOL_SRC_DIR})
  target_link_libraries(${name} PRIVATE ${ARGN})
endfunction()

lillymol_bdb_executable(enough_inventory_bdb
  accumulator iwcmdline iwstring_data_source iw_tdt iwmisc berkeleydb
)

lillymol_bdb_executable(iwbdb_cat
  accumulator iwcmdline iwstring_data_source iwmisc berkeleydb
)

lillymol_bdb_executable(iwbdb_compare
  iwstring_bdb unpack_data accumulator iwcmdline iwstring_data_source iwmisc berkeleydb
)

lillymol_bdb_executable(iwbdb_delete
  accumulator iwcmdline iwstring_data_source iwmisc berkeleydb
)

lillymol_bdb_executable(iwbdb_exists
  accumulator iwcmdline iwstring_data_source iwmisc berkeleydb
)

lillymol_bdb_executable(iwbdb_fetch
  accumulator iwcmdline iwstring_data_source iwbits iwmisc berkeleydb
)

lillymol_bdb_executable(iwbdb_from_tdt
  accumulator iwcmdline iwstring_data_source iw_tdt iwmisc berkeleydb re2::re2 ZLIB::ZLIB
)

lillymol_bdb_executable(iwbdb_list
  unpack_data accumulator iwcmdline iwstring_data_source iwmisc berkeleydb
)

lillymol_bdb_executable(iwbdb_load
  accumulator iwcmdline iwstring_data_source iwmisc berkeleydb
)

lillymol_bdb_executable(iwbdb_merge_into_tdt
  accumulator iwcmdline iwstring_data_source iw_tdt iwmisc berkeleydb re2::re2
)
