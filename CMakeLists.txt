cmake_minimum_required(VERSION 3.14) # version can be different

project(lector_final) #name of your project
set(CMAKE_CXX_STANDARD 14)
add_subdirectory(googletest) # add googletest subdirectory

include_directories(googletest/include) # this is so we can #include <gtest/gtest.h>
include_directories(tests)
include_directories(utils)

add_executable(mytests mytests.cpp utils) # add this executable
add_executable(utils utils/utils.cpp)
add_executable(a Atomo.cpp)
add_executable(main main.cpp)

target_link_libraries(mytests PRIVATE gtest) # link google test to this executable
