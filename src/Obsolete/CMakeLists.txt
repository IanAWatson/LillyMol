add_library(distance_matrix_base
  _dm_to_nn_conditions_float.cc
  _iwdistance_matrix_byte.cc
  _iwdistance_matrix_double.cc
  _iwdistance_matrix_float.cc
  _iwdistance_matrix_float_byte.cc
  _iwdistance_matrix_int.cc
  _masquerading_as_byte_float.cc
  iwdmsupport.cc
)

target_sources(distance_matrix_base
  PUBLIC
    DM_to_NN_Conditions_.h
    IWDistanceMatrixBase.h
    Masquerading_as_Byte.h
    iwdmsupport.h
)

target_include_directories(distance_matrix_base
  PUBLIC
    ${LILLYMOL_SRC_DIR}
)

target_link_libraries(distance_matrix_base
  PUBLIC
    accumulator
    iwstring_data_source
    iwmisc
    iwqsort
)
