project(filters)
set( FILTERS_SOURCES
    minimal_filter_3_moments.cpp
    minimal_filter_factory.cpp
     )
 add_library(filters ${FILTERS_SOURCES})
target_link_libraries(
                        filters
                        utils
                     )
add_clangformat(filters)
add_clangtidy(filters)
