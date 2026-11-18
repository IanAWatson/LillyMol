cmake_minimum_required(VERSION 3.21)
project(LillyMolInchiStatic LANGUAGES C)

set(INCHI_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../InChI" CACHE PATH
    "Path to the upstream InChI checkout")
set(inchi_api "${INCHI_SOURCE_DIR}/INCHI-1-SRC/INCHI_API/libinchi/src")
if(NOT EXISTS "${inchi_api}/CMakeLists.txt")
  message(FATAL_ERROR "InChI sources not found in ${INCHI_SOURCE_DIR}")
endif()

# Upstream explicitly creates a SHARED target, regardless of BUILD_SHARED_LIBS.
# Load its recipe without building that target, then reuse the source list and
# compilation settings for our static archive. No upstream files are modified.
add_subdirectory("${inchi_api}" upstream EXCLUDE_FROM_ALL)
get_target_property(inchi_sources libinchi SOURCES)
get_target_property(inchi_source_dir libinchi SOURCE_DIR)
set(static_sources)
foreach(source IN LISTS inchi_sources)
  if(IS_ABSOLUTE "${source}")
    list(APPEND static_sources "${source}")
  else()
    list(APPEND static_sources "${inchi_source_dir}/${source}")
  endif()
endforeach()

add_library(inchi_static STATIC ${static_sources})
get_target_property(inchi_includes libinchi INCLUDE_DIRECTORIES)
get_target_property(inchi_definitions libinchi COMPILE_DEFINITIONS)
get_target_property(inchi_libraries libinchi LINK_LIBRARIES)
target_include_directories(inchi_static PRIVATE ${inchi_includes})
target_compile_definitions(inchi_static PRIVATE ${inchi_definitions})
target_link_libraries(inchi_static PRIVATE ${inchi_libraries})
set_target_properties(inchi_static PROPERTIES
  OUTPUT_NAME inchi
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
  POSITION_INDEPENDENT_CODE ON)
