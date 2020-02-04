set (NAME Mesh)

set (SOURCES 
Cell.cpp
Face.cpp
Grid.cpp
Line.cpp
Node.cpp
Rectangular.cpp
Structured.cpp
Vertex.cpp
)

add_library(${NAME} OBJECT ${SOURCES} )
 
