add_gmock_test(PointFunctionsTest point_functions_test.cpp)
target_link_libraries(
    PointFunctionsTest
    input_output
    utils
    )
add_clangformat(PointFunctionsTest)

add_gmock_test(FlagHandlerTest flag_handler_test.cpp)
target_link_libraries(
    FlagHandlerTest
    utils
    )
add_clangformat(FlagHandlerTest)

add_gmock_test(ShockDiscontinuityTest shock_discontinuity_test.cpp)
target_link_libraries(
    ShockDiscontinuityTest
    utils
    input_output
    )
add_clangformat(ShockDiscontinuityTest)


