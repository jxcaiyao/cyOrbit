cmake_minimum_required(VERSION 3.10)

# Set the project name
project(cyOrbitTest)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
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

# Add EOP source files
file(GLOB EOP_SOURCES "data/*.cpp")

# Create a static library for SOFA sources
add_library(sofaLib STATIC ${SOFA_SOURCES})
set_target_properties(sofaLib PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

# Add the executables
add_executable(cyOrbitTest src/cyOrbitTest.cpp ${SOURCES} ${LIBS})
add_executable(t_sofa_c src/t_sofa_c.c ${LIBS})

# Add the test executable
add_executable(EOPDataHandlerTest data/EOPDataHandlerTest.cpp ${EOP_SOURCES})

# Link the math library and Google Test
target_link_libraries(cyOrbitTest PRIVATE sofaLib m gtest gtest_main)
target_link_libraries(t_sofa_c PRIVATE sofaLib m)
target_link_libraries(EOPDataHandlerTest PRIVATE sofaLib m gtest gtest_main)

# Set compiler flags based on the build type
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("Build type: Debug")
    target_compile_options(cyOrbitTest PRIVATE -g)
    target_compile_options(t_sofa_c PRIVATE -g)
    target_compile_options(EOPDataHandlerTest PRIVATE -g)
else()
    message("Build type: Release")
    target_compile_options(cyOrbitTest PRIVATE -O2)
    target_compile_options(t_sofa_c PRIVATE -O2)
    target_compile_options(EOPDataHandlerTest PRIVATE -O2)
endif()