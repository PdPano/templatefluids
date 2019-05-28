project(input_output)
set( INPUT_OUTPUT_SOURCES
     options.cpp
     tokenizer.cpp
     stream_from_file.cpp
     )
 add_library(input_output ${INPUT_OUTPUT_SOURCES})
add_subdirectory(readers)
add_subdirectory(writers)
add_subdirectory(test)
target_link_libraries(
                         input_output
                     )
add_clangformat(input_output)
add_clangtidy(input_output)
