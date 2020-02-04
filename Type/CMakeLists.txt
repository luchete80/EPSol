set (NAME Type)

set (SOURCES 
Operations.cpp
Scalar.cpp
Vec3d.cpp
)

add_library(${NAME} OBJECT ${SOURCES} )

install (FILES 
#${NAME}.h 
Matrix.h
Operations.h
Products.h
pTraits.h
Scalar.h
Tensor.h
Vec3d.h
DESTINATION include)
install (FILES lib${NAME}.a DESTINATION lib)
 
