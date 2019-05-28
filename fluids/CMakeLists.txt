cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(templateFluids CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

add_subdirectory(grid)
add_subdirectory(utils)
add_subdirectory(input_output)
add_subdirectory(derivatives)
add_subdirectory(boundary)
add_subdirectory(reconstructions)
add_subdirectory(time_integrators)

add_executable(templatefluids.x main.cpp)
target_link_libraries(
                        templatefluids.x
                        grid
                        input_output
                        time_integrators
                      )
add_clangformat(templatefluids.x)
add_clangtidy(templatefluids.x)
