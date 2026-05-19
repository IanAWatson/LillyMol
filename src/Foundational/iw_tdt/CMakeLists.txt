add_library(iw_tdt
  _iw_tdt_double.cc
  _iw_tdt_float.cc
  _iw_tdt_int.cc
  _iw_tdt_iwstring.cc
  iw_tdt.cc
)

target_include_directories(iw_tdt
  PUBLIC
    ${LILLYMOL_SRC_DIR}
)

target_link_libraries(iw_tdt
  PUBLIC
    iwstring_data_source
    iwaray
    iwmisc
    iwconfig
    iwstring
)

add_executable(tdt_stats tdt_stats.cc)
target_link_libraries(tdt_stats
  PRIVATE
    iw_tdt
    cmdline_v2
    iwmisc
    iwstring
)

add_executable(tdt_sort tdt_sort.cc)
target_link_libraries(tdt_sort
  PRIVATE
    iw_tdt
    cmdline_v2
    iwqsort
    iwmisc
    iwstring
)
