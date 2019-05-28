project(derivatives)
set( DERIVATIVES_SOURCES
    derivatives_factory.cpp
     regular_derivatives.cpp
     regular_derivatives_fourth.cpp
     irregular_derivatives.cpp
     )
add_library(derivatives ${DERIVATIVES_SOURCES})
target_link_libraries(
                         derivatives
                         utils
                         grid
                     )
add_clangformat(derivatives)
add_clangtidy(derivatives)
add_subdirectory(test)
