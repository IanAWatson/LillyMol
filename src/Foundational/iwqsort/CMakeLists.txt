find_package(TBB QUIET)

add_library(iwqsort INTERFACE)

target_sources(iwqsort
  INTERFACE
    iwqsort.h
    iwqsort_ptr.h
    iwqsort_ptr_.h
    iwqsort_tbb.h
)

target_include_directories(iwqsort
  INTERFACE
    ${LILLYMOL_SRC_DIR}
)

target_link_libraries(iwqsort
  INTERFACE
    iwaray
)

if(TBB_FOUND OR TARGET TBB::tbb)
  add_library(iwqsort_parallel INTERFACE)

  target_sources(iwqsort_parallel
    INTERFACE
      iwqsort_tbb.h
  )

  target_include_directories(iwqsort_parallel
    INTERFACE
      ${LILLYMOL_SRC_DIR}
  )

  target_link_libraries(iwqsort_parallel
    INTERFACE
      iwaray
      TBB::tbb
  )
endif()
