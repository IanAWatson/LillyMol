add_library(cmdline_v2
  cmdline.cc
)

target_include_directories(cmdline_v2
  PUBLIC
    ${LILLYMOL_SRC_DIR}
)

target_link_libraries(cmdline_v2
  PUBLIC
    iwaray
    iwmisc
    iwstring
)
