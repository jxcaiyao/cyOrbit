cmake_minimum_required(VERSION 3.10)

# Set the project name
project(cyOrbitTest)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add include directories
include_directories(include lib)

# Add source files
file(GLOB SOURCES "src/*.cpp")

# Add SOFA source files
file(GLOB SOFA_SOURCES "lib/Sofa/*.c")

# Add the executable
add_executable(cyOrbitTest src/cyOrbitTest.cpp ${SOURCES} ${LIBS})
add_executable(t_sofa_c src/t_sofa_c.c ${SOFA_SOURCES} ${LIBS})

target_link_libraries(cyOrbitTest PRIVATE m)
target_link_libraries(t_sofa_c PRIVATE m)