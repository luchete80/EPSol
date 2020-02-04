set (NAME FEGrid)
set (SOURCES 
Element.cpp
FEGrid.cpp
)

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/${NAME}/*.h
  )
  


add_library(${NAME} OBJECT ${SOURCES} )