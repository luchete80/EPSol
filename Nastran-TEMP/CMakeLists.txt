set (NAME Nastran)
set (SOURCES 
Cadenas.cpp
Nastran.cpp
Caract.cpp
Elementos.cpp
Material.cpp
Modelo.cpp
Nastran.cpp
Nodo.cpp
Prop.cpp
SistCoord.cpp
Varios.cpp
)

add_library(${NAME} STATIC ${SOURCES} )

install (FILES ${NAME}.h DESTINATION include)
install (FILES lib${NAME}.a DESTINATION lib)
 