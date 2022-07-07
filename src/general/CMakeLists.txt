#
#  build the files in general folder
#  list all of cpp files
#
set(SRC_GENERAL_CXX 
	${PROJECT_SOURCE_DIR}/src/general/src/basis.cpp  
	${PROJECT_SOURCE_DIR}/src/general/src/basisutil.cpp  
	${PROJECT_SOURCE_DIR}/src/general/src/derivinfor.cpp  
	${PROJECT_SOURCE_DIR}/src/general/src/integral.cpp  
	${PROJECT_SOURCE_DIR}/src/general/src/inttype.cpp  
	${PROJECT_SOURCE_DIR}/src/general/src/kdigestion.cpp  
	${PROJECT_SOURCE_DIR}/src/general/src/jdigestion.cpp  
	${PROJECT_SOURCE_DIR}/src/general/src/shell.cpp  
	${PROJECT_SOURCE_DIR}/src/general/src/shellquartet.cpp
	)
add_library(general ${SRC_GENERAL_CXX})

