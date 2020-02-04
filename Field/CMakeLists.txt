set (NAME Field)

#set (SOURCES 
#Boundary.cpp
#)

#add_library(${NAME} STATIC ${SOURCES} )

install (FILES ${NAME}.h DESTINATION include)
install (FILES lib${NAME}.a DESTINATION lib)
 
