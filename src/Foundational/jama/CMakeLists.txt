add_library(jama INTERFACE)

target_sources(jama
  INTERFACE
    jama_cholesky.h
    jama_eig.h
    jama_lu.h
    jama_qr.h
    jama_svd.h
)

target_include_directories(jama
  INTERFACE
    ${LILLYMOL_SRC_DIR}
)

target_link_libraries(jama
  INTERFACE
    tnt
)
