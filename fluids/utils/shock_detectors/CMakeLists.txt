project(shock_detectors)
set( SHOCK_DETECTORS_SOURCES
    luisa_detector_factory.cpp
    luisa_detector_23.cpp
    luisa_detector_345.cpp
     )
 add_library(shock_detectors ${SHOCK_DETECTORS_SOURCES})
target_link_libraries(
                        shock_detectors
                        utils
                        grid
                        input_output
                        filters
                     )
add_clangformat(shock_detectors)
add_clangtidy(shock_detectors)
add_subdirectory(test)
