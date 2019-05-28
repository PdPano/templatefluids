
project(writers)
set( WRITERS_SOURCES
    writer.cpp
    default_writer.cpp
    vtk_writer.cpp
    nan_checker.cpp
     )
 add_library(writers ${WRITERS_SOURCES})
target_link_libraries(
                        writers
                        input_output
                        grid
                        utils
                     )
add_clangformat(writers)
add_clangtidy(writers)
