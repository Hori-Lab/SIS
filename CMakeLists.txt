cmake_minimum_required(VERSION 3.0)

enable_language(Fortran)
project(sis Fortran)


#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
#endif()

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_subdirectory(src)
