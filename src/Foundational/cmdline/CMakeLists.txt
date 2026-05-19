add_library(iwcmdline
  cmdline.cc
  cmdline_stdstring.cc
)

target_include_directories(iwcmdline
  PUBLIC
    ${LILLYMOL_SRC_DIR}
)

target_link_libraries(iwcmdline
  PUBLIC
    iwconfig
    iwstring
)
