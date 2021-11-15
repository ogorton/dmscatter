cmake_minimum_required(VERSION 3.20.1)
enable_language(Fortran)

project(
    dmfortfactor
    VERSION 1.9
    LANGUAGES Fortran)

option(SERIAL "Option to compile for serial execution." OFF)
if (SERIAL)
    set(CMAKE_Fortran_FLAGS "-Ofast")
else()
    FIND_PACKAGE(OpenMP REQUIRED)
    set(CMAKE_Fortran_FLAGS "-fopenmp -Ofast")
endif()

file(GLOB SOURCES "*.f90")
add_executable(dmfortfactor ${SOURCES})
