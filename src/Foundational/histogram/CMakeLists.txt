add_library(iwhistogram
  histogram_fp.cc
  iwhistogram.cc
)

target_include_directories(iwhistogram
  PUBLIC
    ${LILLYMOL_SRC_DIR}
)

target_link_libraries(iwhistogram
  PUBLIC
    iwbits
    iwstring
)
