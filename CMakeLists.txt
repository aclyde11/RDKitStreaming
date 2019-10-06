cmake_minimum_required( VERSION 3.1 )

set (CMAKE_CXX_STANDARD 14)

project(RDKitSV)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# The compiled library code is here
add_subdirectory(src)

# The executable code is here
add_subdirectory(app)