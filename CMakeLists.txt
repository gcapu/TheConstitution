cmake_minimum_required(VERSION 3.2) 
project(TheConstitution)

option(GOOGLE_TEST "Build google tests" ON)


# The version number. It doesn't do anything yet
set (VERSION_MAJOR 0)
set (VERSION_MINOR 1)     

FIND_PACKAGE( OpenMP )
if(OPENMP_FOUND)
	message("OPENMP FOUND")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if(UNIX)
	include_directories(/opt/local/include/eigen3)
	add_definitions(-std=c++11 -stdlib=libc++ -DNDEBUG -DEIGEN_NO_DEBUG -O3 -g -Wall -DMACOS)
else()
	find_path(EIGEN_INCLUDE_DIR Eigen3)
	if(EIGEN_INCLUDE_DIR STREQUAL "EIGEN_INCLUDE_DIR-NOTFOUND")
	  message(FATAL_ERROR "Could not find Eigen Library.")
	endif()
	include_directories(${EIGEN_INCLUDE_DIR})
endif()


file(GLOB_RECURSE HDRS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h)
set_property( GLOBAL PROPERTY USE_FOLDERS ON)
 
add_library(${PROJECT_NAME} INTERFACE) 
target_sources(${PROJECT_NAME} INTERFACE ${HDRS})
target_include_directories(${PROJECT_NAME} INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:${INCLUDE_INSTALL_DIR}/include/TheConstitution>
  )
 
# tell CMake where to find *.cmake module files
set(CMAKE_MODULE_PATH 
	${CMAKE_MODULE_PATH}
    ${PROJECT_SOURCE_DIR}/cmake
    )

if(GOOGLE_TEST)
	include(tests)
	file(GLOB_RECURSE TEST_SOURCES "tests/*.cpp")
	#add the tests
	add_executable(tests ${TEST_SOURCES})
	target_link_libraries(tests gtest_main TheConstitution)
	#include(CTest)
	enable_testing()
	add_test(NAME dummy_test COMMAND tests)
endif()









 











