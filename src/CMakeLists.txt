cmake_minimum_required(VERSION 3.24)
project(LillyMol LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(LILLYMOL_BUILD_TESTS "Build tests" OFF)
option(LILLYMOL_BUILD_FILECONV "Build fileconv first-pass target" ON)

find_package(ZLIB REQUIRED)
find_package(Protobuf REQUIRED)

find_package(absl QUIET)

set(LILLYMOL_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

set(LILLYMOL_ABSL_CRC32C "")

if(TARGET absl::crc32c)
  set(LILLYMOL_ABSL_CRC32C absl::crc32c)
elseif(TARGET absl::crc)
  set(LILLYMOL_ABSL_CRC32C absl::crc)
else()
  message(STATUS "Abseil CRC target not found; TensorFlow DataRecord support will be skipped")
endif()

if(NOT TARGET absl::status)
  message(FATAL_ERROR "Found Abseil, but target absl::status is missing")
endif()

find_package(tomlplusplus QUIET)

find_path(RE2_INCLUDE_DIR re2/re2.h)
find_library(RE2_LIBRARY NAMES re2)

add_library(re2::re2 UNKNOWN IMPORTED)
set_target_properties(re2::re2 PROPERTIES
  IMPORTED_LOCATION "${RE2_LIBRARY}"
  INTERFACE_INCLUDE_DIRECTORIES "${RE2_INCLUDE_DIR}"
)

option(LILLYMOL_USE_EXISTING_PROTOBUF_FILES
  "Use existing .pb.cc/.pb.h files instead of running protoc"
  ON
)

function(lillymol_proto_library target proto_file)
  get_filename_component(proto_base "${proto_file}" NAME_WE)

  if(LILLYMOL_USE_EXISTING_PROTOBUF_FILES)
    add_library(${target}
      ${CMAKE_CURRENT_SOURCE_DIR}/${proto_base}.pb.cc
      ${CMAKE_CURRENT_SOURCE_DIR}/${proto_base}.pb.h
    )

    target_include_directories(${target}
      PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${LILLYMOL_SRC_DIR}
    )
  else()
    protobuf_generate_cpp(${target}_SRCS ${target}_HDRS ${proto_file})

    add_library(${target}
      ${${target}_SRCS}
      ${${target}_HDRS}
    )

    target_include_directories(${target}
      PUBLIC
        ${CMAKE_CURRENT_BINARY_DIR}
        ${LILLYMOL_SRC_DIR}
    )
  endif()

  target_link_libraries(${target}
    PUBLIC
      ${LILLYMOL_PROTOBUF_LIB}
      ${ARGN}
  )
endfunction()

add_library(lillymol_options INTERFACE)
target_include_directories(lillymol_options INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(lillymol_options INTERFACE ZLIB::ZLIB protobuf::libprotobuf re2::re2)

add_subdirectory(Foundational)
add_subdirectory(Molecule_Lib)

if(LILLYMOL_BUILD_FILECONV)
  add_subdirectory(Molecule_Tools)
endif()


# Not sure if this should be here or on BerkeleyDB
option(LILLYMOL_BUILD_BERKELEYDB "Build Berkeley DB based tools" OFF)

if(LILLYMOL_BUILD_BERKELEYDB)
  find_path(BERKELEYDB_INCLUDE_DIR db.h)
  find_library(BERKELEYDB_LIBRARY NAMES db)

  add_library(berkeleydb UNKNOWN IMPORTED)
  set_target_properties(berkeleydb PROPERTIES
    IMPORTED_LOCATION "${BERKELEYDB_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${BERKELEYDB_INCLUDE_DIR}"
  )
endif()

if(LILLYMOL_BUILD_BERKELEYDB)
  add_subdirectory(BerkeleyDB)
endif()

option(LILLYMOL_BUILD_GFP_KNN "Build legacy GFP_Knn tools" ON)

if(LILLYMOL_BUILD_GFP_KNN)
  add_subdirectory(Utilities/GFP_Knn)
endif()

option(LILLYMOL_FETCH_DEPS "Download missing third-party dependencies" ON)

if(LILLYMOL_BUILD_TESTS)
  enable_testing()
  find_package(GTest QUIET)

  if(NOT TARGET GTest::gtest_main AND LILLYMOL_FETCH_DEPS)
    include(FetchContent)

    FetchContent_Declare(
      googletest
      URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
    )

    # Helps Windows users; harmless elsewhere.
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    FetchContent_MakeAvailable(googletest)
  endif()

  if(NOT TARGET GTest::gtest_main)
    message(FATAL_ERROR
      "GoogleTest not found. Install gtest, configure with "
      "-DLILLYMOL_FETCH_DEPS=ON, or disable tests with "
      "-DLILLYMOL_BUILD_TESTS=OFF"
    )
  endif()
endif()
