set( CARTESIAN_GRID_TEST_INTERFACE_SOURCES
    cartesian_grid_test_interface.cpp
    )
add_library(cartesian_test_interface ${CARTESIAN_GRID_TEST_INTERFACE_SOURCES})
target_link_libraries(
    cartesian_test_interface
    grid
    input_output
    readers
    )
add_clangformat(cartesian_test_interface)

add_gmock_test(CartesianGridTest cartesian_grid_test.cpp)
target_link_libraries(
    CartesianGridTest
    cartesian_test_interface
    grid
    input_output
    readers
    )
add_gmock_test(GhiasGridTest ghias_grid_test.cpp)
target_link_libraries(
    GhiasGridTest
    grid
    input_output
    readers
    utils
    )
add_clangformat(GhiasGridTest)

add_gmock_test(KaragiozisGridTest karagiozis_grid_test.cpp)
target_link_libraries(
    KaragiozisGridTest
    grid
    input_output
    readers
    utils
    )
add_clangformat(KaragiozisGridTest)



add_gmock_test(ShockGridTest shock_grid_test.cpp)
target_link_libraries(
    ShockGridTest
    grid
    input_output
    readers
    utils
    )
add_clangformat(ShockGridTest)
