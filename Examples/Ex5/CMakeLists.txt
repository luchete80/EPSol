cmake_minimum_required (VERSION 2.6)
project (Ex5)

set (CMAKE_CXX_FLAGS "-fpermissive -mwindows")

#include_directories (${SINGLEINPUT_SOURCE_DIR}/Hello) 
include_directories (${EX5_SOURCE_DIR})
link_directories (${EX5_DIR}) 
#/cygdrive/c/cygwin64/petsc-3.4.4
set (PETSC_DIR /petsc-3.4.4)
set (PETSC_CONF arch-mswin-c-release-mpi0-fc0)

#set (PETSC_LIBRARY_DIRS  /media/sf_EPSol_Libs/bin/Ubuntu/petsc)
set (PETSC_LIBRARY_DIRS  ${PETSC_DIR}/${PETSC_CONF}/lib)
#set (PETSC_INCLUDE_DIRS /media/sf_EPSol_Libs/petsc-3.4.4/include)
set (PETSC_INCLUDE_DIRS ${PETSC_DIR}/include)
set (PETSC_INCLUDE_DIRS_MPI ${PETSC_DIR}/include/mpiuni)
set (PETSC_INCLUDE_DIRS_2 ${PETSC_DIR}/${PETSC_CONF}/include)

set (METIS_LIBRARY_DIRS /cygdrive/c/EPSol\ Libs/bin/Win64/Debug/metis-5.1.0-cygwin-mingw-x86-64/libmetis)
set (METIS_INCLUDE_DIRS /cygdrive/c/EPSol\ Libs/metis-5.1.0/include)

include_directories(${METIS_INCLUDE_DIRS})
link_directories(${METIS_LIBRARY_DIRS})


include_directories(${PETSC_INCLUDE_DIRS})
include_directories(${PETSC_INCLUDE_DIRS_MPI})
include_directories(${PETSC_INCLUDE_DIRS_2})
link_directories(${PETSC_LIBRARY_DIRS})


include_directories(/cygdrive/c/Users/samsung/Dropbox/Numerico/MecSol/EPSol/Src)

set (EPSOL_INCLUDE_DIR /cygdrive/c/Users/samsung/Dropbox/Numerico/MecSol/EPSol/Src)

include_directories("${EPSOL_INCLUDE_DIR}/Boundary")
include_directories("${EPSOL_INCLUDE_DIR}/Field")
include_directories("${EPSOL_INCLUDE_DIR}/FiniteElement")
include_directories("${EPSOL_INCLUDE_DIR}/FiniteElement/DoFHandler")
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

set (EPSOL_INCLUDE_DIR /cygdrive/c/Users/samsung/Dropbox/Numerico/MecSol/EPSol/Src)
set (EPSOL_LIB_DIR /cygdrive/c/Users/samsung/Dropbox/Numerico/MecSol/EPSol/Bin/Windows/Debug/EPSol)

link_directories (${EPSOL_LIB_DIR})

set (EXTRA_LIBS petsc f2clapack f2cblas gdi32 metis)
set (EPSOL_LIB epsollib)


#kernel32 user32 gdi32 comdlg32 winspool winmm shell32 comctl32 ole32 oleaut32 uuid rpcrt4 advapi32 wsock32 odbc32 glu32 opengl32

add_executable (Ex5 Ex5.cpp)

#target_link_libraries (SparseExample ${FSOL_LIBS} ${EXTRA_LIBS})

target_link_libraries (Ex5 ${EXTRA_LIBS} ${EPSOL_LIB})
 