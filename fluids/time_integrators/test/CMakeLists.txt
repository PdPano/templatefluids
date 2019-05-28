add_gmock_test(TimeIntegratorTest time_integrator_test.cpp)
target_link_libraries(
    TimeIntegratorTest
    boundary
    derivatives
    cartesian_test_interface
    grid
    input_output
    readers
    utils
    )
add_clangformat(TimeIntegratorTest)
