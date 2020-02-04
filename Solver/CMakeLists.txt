set (NAME Solver)

SET(_src
PETSC_Solver.cpp
)

SET(_inst
PETSC_Solver.inst
)

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/${NAME}/*.h
  )


#OBJECT IS NEEDED TO LINK IN A SINGLE LIBRARY
add_library(${NAME} OBJECT ${_src} ${_inst} ${_header})

#install (FILES ${NAME}.h DESTINATION include)
#install (TARGETS ${NAME} DESTINATION lib)