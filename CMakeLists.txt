cmake_minimum_required(VERSION 3.12)
project(lector-cpmd)
add_subdirectory(tests)
add_subdirectory(utils)

add_executable(angulos angulos.cpp Atomo.cpp Atomo.h utils/utils.cpp utils/utils.h)
add_executable(gdr gdr.cpp Atomo.cpp Atomo.h utils/utils.cpp utils/utils.h)

target_include_directories(angulos PRIVATE utils)
target_include_directories(gdr PRIVATE utils)
target_link_libraries(angulos fftw3)
target_link_libraries(gdr fftw3)
