add_gmock_test(RegularDerivativeTest regular_derivatives_test.cpp)
target_link_libraries(
    RegularDerivativeTest
    derivatives
    cartesian_test_interface
    grid
    input_output
    readers
    utils
    )
add_clangformat(RegularDerivativeTest)

add_gmock_test(IrregularDerivativeTest irregular_derivatives_test.cpp)
target_link_libraries(
    IrregularDerivativeTest
    derivatives
    grid
    input_output
    readers
    utils
    )
add_clangformat(IrregularDerivativeTest)
