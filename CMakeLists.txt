cmake_minimum_required(VERSION 2.8.4)
include (CheckSymbolExists)
include (CheckCXXSymbolExists)
include (CheckIncludeFiles)
include (CheckIncludeFileCXX)
include (CheckCXXSourceCompiles)

project(epsol)

#######################################################################
######################## Options ######################################
#######################################################################
option (USE_PETSC "Enable PETSc support (required for parallel)" ON)
option (USE_METIS "Enable METIS support	                       " ON)
if (${CMAKE_VERSION} VERSION_GREATER "2.8.7") # Compability with older version of cmake
    option (USE_SHARED_LIB "Build shared library" ON)
    option (USE_STATIC_LIB "Build static library" OFF)
endif ()
option (DEBUG_MODE "Enable OOFEG" ON)

if(DEBUG_MODE)
# Default is debug build
#if (NOT CMAKE_BUILD_TYPE)
set (CMAKE_BUILD_TYPE Debug)
else()
	set (CMAKE_BUILD_TYPE Release)
endif ()

#set (CMAKE_CXX_FLAGS "-O3 -fpermissive -g")
set (CMAKE_CXX_FLAGS "-fpermissive -g")
#set (CMAKE_EXE_LINKER_FLAGS "-static-libgcc")

#TO SELECT BUILD TYPE
# CAN BE SET INSTEAD WITH cmake -DCMAKE_BUILD_TYPE=Debug ..
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fpermissive -g -Wall")

#set (LASPACK_LIBRARY_DIRS /home/lucho/PruebaCrossComp/HelloLaspack)


#######################################################################
######################## Custom paths #################################
#######################################################################
#set (MY_PETSC_DIR CACHE PATH "Manual PETSc directory")
#set (MY_PETSC_ARCH CACHE STRING "Manual PETSc arch name (e.g. linux-gnu-cxx-opt)")

set (MY_METIS_DIR D:/Luciano/Numerico/Libs/metis-5.1.0)
set (METIS_LIBRARY_DIRS ${MY_METIS_DIR}/bin-x86_64-mingw)

set (MY_PETSC_DIR "D:/Luciano/Numerico/Libs/petsc-3.4.4")

if(DEBUG_MODE)
set (MY_PETSC_ARCH arch-mswin-c-debug-shared-mpi0-fc0)
else ()
set (MY_PETSC_ARCH arch-mswin-c-release-mpi0-fc0)
endif ()

#######################################################################
######################## External libaries ############################
#######################################################################

## PETSC
if (USE_PETSC)
    find_package (PETSc REQUIRED HINTS "${MY_PETSC_DIR}/${MY_PETSC_ARCH}/conf")
    include_directories ("${MY_PETSC_DIR}/include" "${MY_PETSC_DIR}/${MY_PETSC_ARCH}/include" ${PETSC_INCLUDE_PATH})
	include_directories ("${MY_PETSC_DIR}/include/mpiuni")
    find_library (PETSC_LIB petsc PATH "${MY_PETSC_DIR}/${MY_PETSC_ARCH}/lib/") 
    if (${PETSC_LIB} STREQUAL "PETSC_LIB-NOTFOUND")
        message (FATAL_ERROR "PETSc library not found")
    endif ()
    #set (EXT_LIBS ${EXT_LIBS} ${PETSC_LIB} ${PETSC_PACKAGE_LIBS})
	
	#Blas And Lapack, Manual, can be also automatic
	set (PETSC_LIBRARY_DIRS  ${MY_PETSC_DIR}/${MY_PETSC_ARCH}/lib)
	link_directories(${PETSC_LIBRARY_DIRS})
	set (EXT_LIBS ${EXT_LIBS} ${PETSC_LIB} f2clapack f2cblas)
    
	set (MODULE_LIST ${MODULE_LIST} "PETSc")
    set (USE_MPI OFF)
    add_definitions (-D__PETSC_MODULE)
endif ()

if (USE_METIS)
include_directories ("${MY_METIS_DIR}/include")
link_directories (${METIS_LIBRARY_DIRS})
set (EXT_LIBS ${EXT_LIBS} metis)
endif()

#######################################################################
#######################################################################
#######################################################################

list (APPEND CMAKE_MODULE_PATH "${epsol_SOURCE_DIR}/CMake")

# Temporalmente no lo incluyo
# Pass some CMake settings to source code through a header file
#configure_file (
#	"${PROJECT_SOURCE_DIR}/cmake_vars.h.in"
#	"${PROJECT_BINARY_DIR}/cmake_vars.h"
#)

set (CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/../)

# add to the include search path
include_directories(.)
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/Boundary")
include_directories("${PROJECT_SOURCE_DIR}/Field")
include_directories("${PROJECT_SOURCE_DIR}/FiniteElement")
include_directories("${PROJECT_SOURCE_DIR}/FiniteElement/FEField")
include_directories("${PROJECT_SOURCE_DIR}/FiniteElement/DoFHandler")
include_directories("${PROJECT_SOURCE_DIR}/FiniteElement/Integration")
include_directories("${PROJECT_SOURCE_DIR}/FiniteElement/Mesh")
include_directories("${PROJECT_SOURCE_DIR}/FiniteElement/ShapeFunction")
include_directories("${PROJECT_SOURCE_DIR}/Function")
include_directories("${PROJECT_SOURCE_DIR}/Input")
include_directories("${PROJECT_SOURCE_DIR}/Material")
include_directories("${PROJECT_SOURCE_DIR}/Matrix")
include_directories("${PROJECT_SOURCE_DIR}/Mesh")
include_directories("${PROJECT_SOURCE_DIR}/Model")
include_directories("${PROJECT_SOURCE_DIR}/Output")
include_directories("${PROJECT_SOURCE_DIR}/SistEcuac")
include_directories("${PROJECT_SOURCE_DIR}/Solver")
include_directories("${PROJECT_SOURCE_DIR}/Type")
include_directories("${PROJECT_SOURCE_DIR}/Utils")


add_subdirectory(Boundary)
add_subdirectory(Field)
#add_subdirectory(FiniteElement)
add_subdirectory(FiniteElement/Mesh)
add_subdirectory(Input)
add_subdirectory(Mesh)
add_subdirectory(Type)
add_subdirectory(Utils)
add_subdirectory(Solver)
#Pending for output
#set (FSOL_LIBS Boundary Mesh Input Type Utils)

#set (EXTRA_LIBS Laspack)

#add the executable
set (SOURCES
EPSol.cpp
)



set (NAME epsol)
#install (TARGETS epsol RUNTIME DESTINATION bin)

#EPSOL LIBRARY
#add_library(${NAME} STATIC ${SOURCES})
# NO EXTRA SOURCES

if (${CMAKE_VERSION} VERSION_LESS "2.8.8")
    set (USE_SHARED_LIB OFF)
endif ()
if (USE_SHARED_LIB)
#Before creating library
set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--export-all-symbols")

add_library(${NAME}lib SHARED EPSol_lib.cpp 	$<TARGET_OBJECTS:Boundary>
											$<TARGET_OBJECTS:FEGrid>
											$<TARGET_OBJECTS:Mesh>
											$<TARGET_OBJECTS:Input> 
											$<TARGET_OBJECTS:Type>
											$<TARGET_OBJECTS:Solver>
											$<TARGET_OBJECTS:Utils>)
#target_link_libraries (${NAME} ${EXTRA_LIBS} ${EXT_LIBS})
target_link_libraries (${NAME}lib ${EXTRA_LIBS} ${EXT_LIBS})
else ()
    set (USE_STATIC_LIB ON)
endif ()

if (USE_STATIC_LIB)
add_library(${NAME}lib_static STATIC EPSol_lib.cpp 	$<TARGET_OBJECTS:Boundary>
											$<TARGET_OBJECTS:FEGrid>
											$<TARGET_OBJECTS:Mesh>
											$<TARGET_OBJECTS:Input> 
											$<TARGET_OBJECTS:Type>
											$<TARGET_OBJECTS:Solver>
											$<TARGET_OBJECTS:Utils>)
#target_link_libraries (${NAME} ${EXTRA_LIBS} ${EXT_LIBS})
target_link_libraries (${NAME}lib_static ${EXTRA_LIBS} ${EXT_LIBS})
endif()

#add_library(${NAME} SHARED $<TARGET_OBJECTS:subproject1>
#                                   $<TARGET_OBJECTS:subproject2>)


#set_target_properties (${NAME} PROPERTIES OUTPUT_NAME $NAME)

install (TARGETS ${NAME}lib DESTINATION lib)
install (FILES ${NAME}.h DESTINATION include)

#EPSOL EXECUTABLE
add_executable(${NAME} ${SOURCES})
target_link_libraries (${NAME} ${NAME}lib)

#TO MODIFY
#if (USE_SHARED_LIB) # Prefer dynamic lib if available.
#    target_link_libraries (oofem liboofem)
#target_link_libraries (${NAME} ${FSOL_LIBS} ${EXTRA_LIBS})
#install (TARGETS ${NAME} DESTINATION RUNTIME)


#install (FILES ${NAME}.h DESTINATION include)
#install (FILES lib${NAME}.a DESTINATION lib)
