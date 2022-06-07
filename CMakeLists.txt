################################################
#       heading for the top cmake file         #
################################################
cmake_minimum_required(VERSION 3.12)
project(cppints LANGUAGES CXX)
set(BIN_NAME "cppints.exe")

################################################
#        set the compiler flags                #
################################################
set(CMAKE_CXX_FLAGS_RELEASE " ")
set(CMAKE_CXX_FLAGS_DEBUG   " ")

# debug option
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")

# release options
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11")


################################################
#            set the library                   #
# here inside we only use boost                #
################################################
set(THIS_PROJECT_LIBS "")

#
# boost library
#
set(Boost_USE_STATIC_LIBS     OFF)
set(Boost_USE_MULTITHREADED   OFF)
set(Boost_USE_STATIC_RUNTIME  OFF)
find_package(Boost 1.60.0 COMPONENTS filesystem system)
if(Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
	set(THIS_PROJECT_LIBS ${Boost_LIBRARIES} ${THIS_PROJECT_LIBS})
	
	# output message
	message("boost libs: " ${Boost_LIBRARIES})
	message("boost include dir: " ${Boost_INCLUDE_DIRS})
	message("          ")
else()
	message(FATAL_ERROR "boost library not found")
endif()


################################################
#           set the include dir                #
################################################

include_directories ("${PROJECT_SOURCE_DIR}/src/general/include")

################################################
#           add in each subdir                 #
################################################

add_subdirectory("${PROJECT_SOURCE_DIR}/src/general")
set (MY_WHOLE_LIBRARIES "general")

################################################
#     now we need to add the target module     #
################################################

# add exe
add_executable(${BIN_NAME} ${PROJECT_SOURCE_DIR}/src/main.cpp)

# set the linking property
# we use the C++ as linker
# also set the linker option according to the compiler
set_target_properties(${BIN_NAME} PROPERTIES LINKER_LANGUAGE CXX)

# form the target linking library
# THIS_PROJECT_LIBS refers to the third party libraries like boost, mkl used
# in this code
target_link_libraries(${BIN_NAME} ${MY_WHOLE_LIBRARIES} ${THIS_PROJECT_LIBS})
