cmake_minimum_required (VERSION 2.6)
project (Ex3)

set (CMAKE_CXX_FLAGS "-fpermissive -mwindows")

#include_directories (${SINGLEINPUT_SOURCE_DIR}/Hello) 
include_directories (${EX3_SOURCE_DIR})
link_directories (${EX3_DIR}) 
#/cygdrive/c/cygwin64/petsc-3.4.4
set (PETSC_DIR /petsc-3.4.4)
set (PETSC_CONF arch-mswin-c-debug-shared-mpi0-fc0)

#set (PETSC_LIBRARY_DIRS  /media/sf_EPSol_Libs/bin/Ubuntu/petsc)
set (PETSC_LIBRARY_DIRS  ${PETSC_DIR}/${PETSC_CONF}/lib)
#set (PETSC_INCLUDE_DIRS /media/sf_EPSol_Libs/petsc-3.4.4/include)
set (PETSC_INCLUDE_DIRS ${PETSC_DIR}/include)
set (PETSC_INCLUDE_DIRS_MPI ${PETSC_DIR}/include/mpiuni)
set (PETSC_INCLUDE_DIRS_2 ${PETSC_DIR}/${PETSC_CONF}/include)

link_directories(${PETSC_LIBRARY_DIRS})

include_directories(${PETSC_INCLUDE_DIRS})
include_directories(${PETSC_INCLUDE_DIRS_MPI})
include_directories(${PETSC_INCLUDE_DIRS_2})

set (EXTRA_LIBS petsc f2clapack f2cblas gdi32)


#kernel32 user32 gdi32 comdlg32 winspool winmm shell32 comctl32 ole32 oleaut32 uuid rpcrt4 advapi32 wsock32 odbc32 glu32 opengl32

add_executable (Ex3 Ex3.cpp)

#target_link_libraries (SparseExample ${FSOL_LIBS} ${EXTRA_LIBS})

target_link_libraries (Ex3 ${EXTRA_LIBS})
 