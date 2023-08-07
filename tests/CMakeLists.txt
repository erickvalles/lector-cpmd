cmake_minimum_required(VERSION 3.12)
project(tests)

# Agrega los archivos fuente para los tests
add_executable(tests test_calculaPosicionesPeriodicas.cpp ../Atomo.cpp)

# Agrega el enlace a la biblioteca creada en el subproyecto "utils"
target_link_libraries(tests utils)

# Agrega Catch2 si lo has instalado y configurado
include_directories(/ruta/a/la/carpeta/catch2)
find_package(Catch2 REQUIRED)
target_link_libraries(tests Catch2::Catch2)

# Indica el directorio de salida para los archivos binarios de tests
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests)
target_link_libraries(tests fftw3)