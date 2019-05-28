find_program(CLANGFORMAT_EXEC clang-format)
if(CLANGFORMAT_EXEC)
    file(GLOB_RECURSE ALL_SOURCE_FILES *.cpp *.hpp)

    add_custom_target(
        clangformat
        COMMAND ${CLANGFORMAT_EXEC}
        -style=file
        -i
        ${ALL_SOURCE_FILES}
        )
endif()
