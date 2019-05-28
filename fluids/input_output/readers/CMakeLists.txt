project(readers)
set( READERS_SOURCES
    reader.cpp
    default_reader.cpp
     )
 add_library(readers ${READERS_SOURCES})
target_link_libraries(
                        readers
                        input_output
                        utils
                     )
add_clangformat(readers)
add_clangtidy(readers)
