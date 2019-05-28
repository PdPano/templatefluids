SET(CMAKE_CXX_COMPILER "/usr/bin/g++" CACHE string "g++ compiler" FORCE)
SET( CMAKE_CXX_FLAGS_OPT2 "-O2 -fopenmp" CACHE STRING
    "Flags used by the C++ compiler during opt2 builds."
    FORCE )
SET( CMAKE_C_FLAGS_OPT2 "-O2 -fopenmp" CACHE STRING
    "Flags used by the C compiler during opt2 builds."
    FORCE )
SET( CMAKE_CXX_FLAGS_OPT3 "-O3 -fopenmp -g" CACHE STRING
    "Flags used by the C++ compiler during opt3 builds."
    FORCE )
SET( CMAKE_C_FLAGS_OPT3 "-O3 -fopenmp -g" CACHE STRING
    "Flags used by the C compiler during opt3 builds."
    FORCE )
# Update the documentation string of CMAKE_BUILD_TYPE for GUIs
 SET( CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
     "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel Opt2."
         FORCE )
