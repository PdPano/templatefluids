add_gmock_test(FluxFunctionsTest flux_functions_test.cpp)
target_link_libraries(
    FluxFunctionsTest
    input_output
    utils
    )
add_clangformat(FluxFunctionsTest)

