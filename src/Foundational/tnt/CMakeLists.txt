add_library(tnt INTERFACE)

target_sources(tnt
  INTERFACE
    tnt_array1d.h
    tnt_array2d.h
    tnt_array3d.h
    tnt_cmat.h
    tnt_fortran_array1d.h
    tnt_fortran_array2d.h
    tnt_fortran_array3d.h
    tnt_i_refvec.h
    tnt_math_utils.h
    tnt_sparse_matrix_csr.h
    tnt_stopwatch.h
    tnt_subscript.h
    tnt_vec.h
)

target_include_directories(tnt
  INTERFACE
    ${LILLYMOL_SRC_DIR}
)
