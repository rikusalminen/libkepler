cmake_minimum_required(VERSION 2.8)
project(libkepler C)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type for make-style systems" FORCE)
endif()

if (CMAKE_COMPILER_IS_GNUCC)
    add_definitions(-std=gnu99)
    add_definitions(-W -Wall -Wextra)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_subdirectory(src)
add_subdirectory(test)
