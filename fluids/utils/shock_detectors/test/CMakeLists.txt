add_gmock_test(LuisaDetectorTest luisa_detector_test.cpp)
target_link_libraries(
    LuisaDetectorTest
    utils
    input_output
    cartesian_test_interface
    )
