cmake_minimum_required(VERSION 2.8)
include (CheckSymbolExists)
include (CheckCXXSymbolExists)
include (CheckIncludeFiles)
include (CheckIncludeFileCXX)
include (CheckCXXSourceCompiles)

project(Ex1)

#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -fpermissive")
set (CMAKE_CXX_FLAGS "-fpermissive -static-libgcc -static-libstdc++")


# add to the include search path
include_directories(.)
#include_directories("${PROJECT_SOURCE_DIR}/../..")
include_directories(/cygdrive/c/Users/samsung/Dropbox/Numerico/MecSol/EPSol/Src)

set (EPSOL_INCLUDE_DIR /cygdrive/c/Users/samsung/Dropbox/Numerico/MecSol/EPSol/Src)
set (EPSOL_LIB_DIR /cygdrive/c/Users/samsung/Dropbox/Numerico/MecSol/EPSol/Bin/Windows/Debug/EPSol)

set (LASPACK_INCLUDE /cygdrive/c/EPSol\ Libs/Laspack/src/)
set (LASPACK_LIBRARY_DIR /cygdrive/c/EPSol\ Libs/Laspack/bin/Win64)

include_directories("${EPSOL_INCLUDE_DIR}/Boundary")
include_directories("${EPSOL_INCLUDE_DIR}/Field")
include_directories("${EPSOL_INCLUDE_DIR}/FiniteElement")
include_directories("${EPSOL_INCLUDE_DIR}/FiniteElement/FEField")
include_directories("${EPSOL_INCLUDE_DIR}/FiniteElement/Integration")
include_directories("${EPSOL_INCLUDE_DIR}/FiniteElement/Mesh")
include_directories("${EPSOL_INCLUDE_DIR}/FiniteElement/ShapeFunction")
include_directories("${EPSOL_INCLUDE_DIR}/Input")
include_directories("${EPSOL_INCLUDE_DIR}/Matrix")
include_directories("${EPSOL_INCLUDE_DIR}/Mesh")
include_directories("${EPSOL_INCLUDE_DIR}/Model")
include_directories("${EPSOL_INCLUDE_DIR}/Output")
include_directories("${EPSOL_INCLUDE_DIR}/SistEcuac")
include_directories("${EPSOL_INCLUDE_DIR}/Solver")
include_directories("${EPSOL_INCLUDE_DIR}/Type")
include_directories("${EPSOL_INCLUDE_DIR}/Utils")

include_directories(/cygdrive/c/EPSol\ Libs/Laspack/src/)
include_directories(${LASPACK_INCLUDE})
link_directories(${LASPACK_LIBRARY_DIR})
link_directories (${EPSOL_LIB_DIR})



set (EPSOL_LIB epsol)

#add the executable
set (SOURCES
Ex4.cpp SingleInputFile.cpp
)
add_executable(Ex1 ${SOURCES})

target_link_libraries (Ex1 ${EPSOL_LIB} Laspack)
