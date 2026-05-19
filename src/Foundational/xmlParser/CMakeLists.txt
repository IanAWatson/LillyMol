add_library(xmlparser
  xmlParser.cpp
)

target_include_directories(xmlparser
  PUBLIC
    ${LILLYMOL_SRC_DIR}
)
