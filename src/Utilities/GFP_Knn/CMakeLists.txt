add_library(extract_from_tdt_form
  extract_from_tdt_form.cc
)

target_include_directories(extract_from_tdt_form PUBLIC ${LILLYMOL_SRC_DIR})

target_link_libraries(extract_from_tdt_form
  PUBLIC
    iwstring
)

add_library(nn_support
  nn_specification.cc
  nn_stuff.cc
  prediction_weight_function.cc
)

target_include_directories(nn_support PUBLIC ${LILLYMOL_SRC_DIR})

target_link_libraries(nn_support
  PUBLIC
    extract_from_tdt_form
    iwbits
    iwmisc
    distance_matrix_base
    gfp
)

add_library(nn_specification
  nn_specification.cc
)

target_include_directories(nn_specification PUBLIC ${LILLYMOL_SRC_DIR})

target_link_libraries(nn_specification
  PUBLIC
    nn_support
    accumulator
    iwcmdline
    iwstring_data_source
    iw_tdt
    iwmisc
)

add_executable(nn_predictions
  nn_predictions.cc
)

target_link_libraries(nn_predictions
  PRIVATE
    nn_specification
    nn_support
    accumulator
    iwcmdline
    iwstring_data_source
    iw_tdt
    iwmisc
    gfp
)

add_executable(nn_training
  nn_training.cc
)

target_link_libraries(nn_training
  PRIVATE
    extract_from_tdt_form
    iwpvalue
    nn_support
    accumulator
    iwcmdline
    iwstring_data_source
    iw_tdt
    iwmisc
    distance_matrix_base
    bsquared
    gfp
)
