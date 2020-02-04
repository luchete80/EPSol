set (NAME Analysis)

set (SOURCES 
Analysis.cpp
)

add_library(${NAME} STATIC ${SOURCES} )

install (FILES ${NAME}.h DESTINATION include)
install (FILES lib${NAME}.a DESTINATION lib)
 
