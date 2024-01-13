project(
  "mt_stream"
  LANGUAGES "Fortran"
  VERSION "1.11"
  DESCRIPTION "Multiple stream Mersenne Twister PRNG"
)

option(MT_USE_NORANGECHECK "Use -fno-range-check in MT module" OFF)
option(MT_USE_ALLOWINVALIDBOZ "Use -fallow-invalid-boz in MT module" ON)

if (MT_USE_NORANGECHECK)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-range-check")
endif()
if (MT_USE_ALLOWINVALIDBOZ)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-invalid-boz")
endif()

add_library("mt_stream-lib" STATIC mt_stream.F90 mt_kind_defs.F90 f_get_coeff.F90 gf2xe.F90 f_jump_coeff.F90)

#add_library("${PROJECT_NAME}" INTERFACE)
#target_link_libraries("${PROJECT_NAME}" INTERFACE "${PROJECT_NAME}-lib")

set_target_properties(
  "mt_stream-lib"
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME "${PROJECT_NAME}"
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
)

target_include_directories(
  "mt_stream-lib"
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}>
)
