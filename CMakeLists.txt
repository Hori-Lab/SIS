cmake_minimum_required(VERSION 3.0)

enable_language(Fortran)
project(sis Fortran)


#ifort
if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    #set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -traceback -check all -warn all -fp-model strict -standf18")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -traceback -check all -warn all -fp-model strict")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O3 -g -march=native")
endif()

#gfortran
if (CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -Wall -g -fbacktrace -ffpe-trap=zero,overflow,denormal,underflow -ffree-line-length-none -fbounds-check -fmax-errors=5")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -ffree-line-length-none")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O3 -g march=native -ffree-line-length-none")
endif()


find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()


set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


# To output Git commit hash. Either CheckGitSetup or NoGitSetup will be called.
include(./CheckGit.cmake)

find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")

    # Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()

    # Output Git hash
    CheckGitSetup()

else()

    # Set githash = "?" and output version number (in main.F90)
    NoGitSetup()

endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/subprojects/toml-f/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()

add_subdirectory(subprojects)
add_subdirectory(src)
