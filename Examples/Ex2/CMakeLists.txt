cmake_minimum_required(VERSION 2.8)
include (CheckSymbolExists)
include (CheckCXXSymbolExists)
include (CheckIncludeFiles)
include (CheckIncludeFileCXX)
include (CheckCXXSourceCompiles)

project(Ex1)

#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -fpermissive")
set (CMAKE_CXX_FLAGS "-fpermissive")


# add to the include search path
include_directories(.)

set (EPSOL_INCLUDE_DIR "$ENV{EPSOL_INCLUDE_DIR}")
set (EPSOL_LIB_DIR "$ENV{EPSOL_LIB_DIR}")

set (LASPACK_INCLUDE /cygdrive/c/EPSol\ Libs/Laspack/src/)
set (LASPACK_LIBRARY_DIR /cygdrive/c/EPSol\ Libs/Laspack/bin/Win64)
set (MY_METIS_DIR "$ENV{METIS_DIR}")
set (METIS_LIBRARY_DIRS "${MY_METIS_DIR}/bin-x86_64-mingw")
set (PETSC_DIR "$ENV{PETSC_DIR}")
set (PETSC_ARCH "$ENV{PETSC_ARCH}")
set (MPI_LIB_DIR "$ENV{MPI_LIB_DIR}")
set (MPI_INC_DIR "$ENV{MPI_INC_DIR}")

include_directories("${EPSOL_INCLUDE_DIR}")
include_directories("${EPSOL_INCLUDE_DIR}/Boundary")
include_directories("${EPSOL_INCLUDE_DIR}/Field")
include_directories("${EPSOL_INCLUDE_DIR}/FiniteElement")
include_directories("${EPSOL_INCLUDE_DIR}/FiniteElement/DoFHandler")
include_directories("${EPSOL_INCLUDE_DIR}/FiniteElement/FEField")
include_directories("${EPSOL_INCLUDE_DIR}/FiniteElement/Integration")
include_directories("${EPSOL_INCLUDE_DIR}/FiniteElement/Mesh")
include_directories("${EPSOL_INCLUDE_DIR}/FiniteElement/ShapeFunction")
include_directories("${EPSOL_INCLUDE_DIR}/Function")
include_directories("${EPSOL_INCLUDE_DIR}/Input")
include_directories("${EPSOL_INCLUDE_DIR}/Material")
include_directories("${EPSOL_INCLUDE_DIR}/Matrix")
include_directories("${EPSOL_INCLUDE_DIR}/Mesh")
include_directories("${EPSOL_INCLUDE_DIR}/Model")
include_directories("${EPSOL_INCLUDE_DIR}/Output")
include_directories("${EPSOL_INCLUDE_DIR}/SistEcuac")
include_directories("${EPSOL_INCLUDE_DIR}/Solver")
include_directories("${EPSOL_INCLUDE_DIR}/Type")
include_directories("${EPSOL_INCLUDE_DIR}/Utils")

include_directories(/cygdrive/c/EPSol\ Libs/Laspack/src/)
include_directories ("${MY_METIS_DIR}/include")
include_directories(${LASPACK_INCLUDE})
link_directories(${LASPACK_LIBRARY_DIR})
link_directories (${EPSOL_LIB_DIR})
link_directories ("${METIS_LIBRARY_DIRS}")
include_directories ("${PETSC_DIR}/include")
include_directories ("${PETSC_DIR}/${PETSC_ARCH}/include")
include_directories ("${PETSC_DIR}/include/mpiuni")
include_directories(${MPI_INC_DIR}

set (EXTRA_LIBS petsc flapack fblas metis msmpi gfortran))

set (EPSOL_LIB epsollib metis)

#add the executable
set (SOURCES
Ex2.cpp
)
add_executable(Ex2 ${SOURCES})

target_link_libraries (Ex2 ${EPSOL_LIB} ${EXTRA_LIBS})
