project(flux_functions)
set( FLUX_FUNCTIONS_SOURCES
    flux_factory.cpp
    simple_flux.cpp
    steger_warming_flux.cpp
    LF_flux.cpp
     )
 add_library(flux_functions ${FLUX_FUNCTIONS_SOURCES})
target_link_libraries(
                        flux_functions
                        utils
                        input_output
                     )
add_clangformat(flux_functions)
add_clangtidy(flux_functions)
