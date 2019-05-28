project(utils)
set( UTILS_SOURCES
     point_functions.cpp
     useful_alias.cpp
     flag_handler.cpp
     operators_overloads.cpp
     shock_discontinuity_handler.cpp
     global_vars.cpp
     )
 add_library(utils ${UTILS_SOURCES})
add_clangformat(utils)
add_clangtidy(utils)
add_subdirectory(shock_detectors)
add_subdirectory(filters)
add_subdirectory(test)
