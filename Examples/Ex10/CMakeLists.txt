cmake_minimum_required (VERSION 2.6)
project (Ex10)

set (CMAKE_BUILD_TYPE Release)

#set (CMAKE_CXX_FLAGS "-fpermissive -mwindows")
set (CMAKE_CXX_FLAGS "-fpermissive")

#include_directories (${SINGLEINPUT_SOURCE_DIR}/Hello) 
include_directories (${Ex10_SOURCE_DIR})
link_directories (${Ex10_DIR}) 


set (PETSC_DIR "D:/Luciano/Numerico/Libs/petsc-3.4.4")
set (PETSC_CONF arch-mswin-c-release-mpi0-fc0)

set (PETSC_LIBRARY_DIRS  ${PETSC_DIR}/${PETSC_CONF}/lib)
set (PETSC_INCLUDE_DIRS ${PETSC_DIR}/include)
set (PETSC_INCLUDE_DIRS_MPI ${PETSC_DIR}/include/mpiuni)
set (PETSC_INCLUDE_DIRS_2 ${PETSC_DIR}/${PETSC_CONF}/include)

set (MY_METIS_DIR D:/Luciano/Numerico/Libs/metis-5.1.0)
set (METIS_LIBRARY_DIRS ${MY_METIS_DIR}/bin-x86_64-mingw)
set (METIS_INCLUDE_DIRS ${MY_METIS_DIR}/include)

link_directories(${PETSC_LIBRARY_DIRS})
link_directories(${METIS_LIBRARY_DIRS})

include_directories(${PETSC_INCLUDE_DIRS})
include_directories(${PETSC_INCLUDE_DIRS_MPI})
include_directories(${PETSC_INCLUDE_DIRS_2})
include_directories(${METIS_INCLUDE_DIRS})

#only for mg2metis
include_directories(${MY_METIS_DIR}/programs)
include_directories(${MY_METIS_DIR}/GKlib)

set (EPSOL_DIR D:/Luciano/Numerico/MecSol/EPSol)
set (EPSOL_INCLUDE_DIR ${EPSOL_DIR}/Src)
set (EPSOL_LIB_DIR ${EPSOL_DIR}/Bin_CodeBlocks/EPSol)

include_directories(${EPSOL_INCLUDE_DIR})
link_directories (${EPSOL_LIB_DIR})


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
include_directories("${EPSOL_INCLUDE_DIR}/Material")
include_directories("${EPSOL_INCLUDE_DIR}/Mesh")
include_directories("${EPSOL_INCLUDE_DIR}/Model")
include_directories("${EPSOL_INCLUDE_DIR}/Output")
include_directories("${EPSOL_INCLUDE_DIR}/SistEcuac")
include_directories("${EPSOL_INCLUDE_DIR}/Solver")
include_directories("${EPSOL_INCLUDE_DIR}/Type")
include_directories("${EPSOL_INCLUDE_DIR}/Utils")


set (EXTRA_LIBS petsc f2clapack f2cblas gdi32 metis)
set (EPSOL_LIB epsollib)

#kernel32 user32 gdi32 comdlg32 winspool winmm shell32 comctl32 ole32 oleaut32 uuid rpcrt4 advapi32 wsock32 odbc32 glu32 opengl32

add_executable (Ex10 Ex10.cpp)

#target_link_libraries (SparseExample ${FSOL_LIBS} ${EXTRA_LIBS})

target_link_libraries (Ex10 ${EPSOL_LIB} ${EXTRA_LIBS})
set_target_properties(Ex10 PROPERTIES PREFIX "") 
install(TARGETS Ex10		RUNTIME DESTINATION bin COMPONENT libraries 
						LIBRARY DESTINATION lib COMPONENT libraries 
						#ARCHIVE DESTINATION lib/static COMPONENT libraries) 
						)
 