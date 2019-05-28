project(time_integrators)
set( TIME_INTEGRATORS_SOURCES
     time_integrator_tool.cpp
     time_integrator_tool_factory.cpp
     solver.cpp
     )
 add_library(time_integrators ${TIME_INTEGRATORS_SOURCES})
target_link_libraries(
                         time_integrators
                         boundary
                         utils
                         grid
                         derivatives
                         convection
                         input_output
                         writers
                         filters
                     )
add_clangformat(time_integrators)
add_clangtidy(time_integrators)
#add_subdirectory(test)
