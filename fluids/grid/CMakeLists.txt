project(grid)
set( GRID_SOURCES
     cartesian_grid.cpp
     ghias_grid.cpp
     ghias_shock_grid.cpp
     karagiozis_grid.cpp
     shock_grid.cpp
     )
add_library(grid ${GRID_SOURCES})
target_link_libraries(
                         grid
                         input_output
                         readers
                         lapack
                         blas
                         shock_detectors
                         utils
                     )
add_clangformat(grid)
add_clangtidy(grid)
add_subdirectory(test)
