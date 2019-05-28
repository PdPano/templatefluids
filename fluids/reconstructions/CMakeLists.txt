project(convection)
add_subdirectory(flux_functions)
set( CONVECTION_SOURCES
     skew_symmetric.cpp
     simple_convection.cpp
     dissipation_tool.cpp
     simple_dissipation.cpp
     split_convection.cpp
     split_convection_cached.cpp
     zero_dissipation.cpp
     dissipation_factory.cpp
     convection_factory.cpp
     simple_flux_convection.cpp
     weno_convection.cpp
     mix_convection.cpp
     )
 add_library(convection ${CONVECTION_SOURCES})
target_link_libraries(
                         convection
                         derivatives
                         utils
                         grid
                         flux_functions
                     )
add_clangformat(convection)
add_clangtidy(convection)
#add_subdirectory(test)
