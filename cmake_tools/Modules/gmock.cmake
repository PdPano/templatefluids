#Based on source from https://www.johnlamp.net/cmake-tutorial.html

set(GMOCK_DIR "./googletest"
    CACHE PATH "The path to the GoogleMock test framework.")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # force this option to ON so that Google Test will use /MD instead of /MT
    # /MD is now the default for Visual Studio, so it should be our default, too
    option(gtest_force_shared_crt
           "Use shared (DLL) run-time lib even when Google Test is built as static lib."
           ON)
elseif (APPLE)
    add_definitions(-DGTEST_USE_OWN_TR1_TUPLE=1)
endif()
add_subdirectory(${GMOCK_DIR} ${CMAKE_BINARY_DIR}/gmock)
set_property(TARGET gtest APPEND_STRING PROPERTY COMPILE_FLAGS " -w")

include_directories(SYSTEM ${GMOCK_DIR}/googletest/include
                           ${GMOCK_DIR}/include/gmock)


#
# add_gmock_test(<target> <sources>...)
#
#  Adds a Google Mock based test executable, <target>, built from <sources> and
#  adds the test so that CTest will run it. Both the executable and the test
#  will be named <target>.
#
function(add_gmock_test_base target source)
    add_executable(${target} ${source})

    add_test(${target} ${target}) 

    add_custom_command(TARGET ${target}
                       POST_BUILD
                       COMMAND ${target} ${ARGN}
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       COMMENT "Running ${target}" VERBATIM)
endfunction()

function(add_gmock_main_test target source)
    add_gmock_test_base(${target} ${source} ${ARGN})
    target_link_libraries(${target} gmock)
endfunction()

function(add_gmock_test target source)
    add_gmock_test_base(${target} ${source} ${ARGN})
    target_link_libraries(${target} gmock_main)
endfunction()
