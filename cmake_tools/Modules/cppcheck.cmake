find_program(CPPCHECK_EXEC cppcheck)
if(CPPCHECK_EXEC)
    file(GLOB_RECURSE ALL_SOURCE_FILES *.cpp *.hpp)
    add_custom_target(
        cppcheck
        COMMAND ${CPPCHECK_EXEC}
        --enable=warning,performance,portability,missingInclude
        --std=c++11
        --verbose
        --quiet
        --check-config
        ${ALL_SOURCE_FILES}
        )
endif()
