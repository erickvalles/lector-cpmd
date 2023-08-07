cmake_minimum_required(VERSION 3.12)
project(utils)

# Agrega los archivos fuente para el subproyecto "utils"
add_library(utils utils.cpp)

# Si tienes otras bibliotecas o dependencias, agrega sus enlaces aquí

# Indica el directorio de salida para los archivos binarios de la biblioteca
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
