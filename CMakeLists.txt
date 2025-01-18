cmake_minimum_required(VERSION 3.10)

# Set the project name
project(cyOrbitTest)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set the default build type to Debug if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# Add include directories
include_directories(include lib)

# Add source files
file(GLOB SOURCES "src/*.cpp")

# Add SOFA source files
file(GLOB SOFA_SOURCES "lib/Sofa/*.c")

# Add the executables
add_executable(cyOrbitTest src/cyOrbitTest.cpp ${SOFA_SOURCES} ${SOURCES} ${LIBS})
add_executable(t_sofa_c src/t_sofa_c.c ${SOFA_SOURCES} ${LIBS})

# Link the math library
target_link_libraries(cyOrbitTest PRIVATE m)
target_link_libraries(t_sofa_c PRIVATE m)

# Set compiler flags based on the build type
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("Build type: Debug")
    target_compile_options(cyOrbitTest PRIVATE -g)
    target_compile_options(t_sofa_c PRIVATE -g)
else()
    message("Build type: Release")
    target_compile_options(cyOrbitTest PRIVATE -O3)
    target_compile_options(t_sofa_c PRIVATE -O3)
endif()