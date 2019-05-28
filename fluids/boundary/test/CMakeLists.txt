add_gmock_test(BoundaryTest boundary_test.cpp)
target_link_libraries(
    BoundaryTest
    boundary
    derivatives
    cartesian_test_interface
    grid
    input_output
    readers
    utils
    )
add_clangformat(BoundaryTest)
