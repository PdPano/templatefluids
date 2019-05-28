project(boundary)
set( BOUNDARY_SOURCES
     boundary.cpp
     lodi.cpp
     )
 add_library(boundary ${BOUNDARY_SOURCES})
target_link_libraries(
                         boundary
                         utils
                         grid
                         convection
                     )
add_clangformat(boundary)
add_clangtidy(boundary)
add_subdirectory(test)
