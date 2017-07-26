
find_package(GTest)
if(${GTEST_FOUND})
	include_directories(${GTEST_INCLUDE_DIRS})
else()
	message("Downloading and building google test.")
	# Download and unpack googletest at configure time
	configure_file(cmake/CMakeLists.txt.in googletest-download/CMakeLists.txt)
	execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
	  RESULT_VARIABLE result
	  OUTPUT_QUIET
	  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
	if(result)
	  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
	endif()
	execute_process(COMMAND ${CMAKE_COMMAND} --build .
	  RESULT_VARIABLE result
	  OUTPUT_QUIET
	  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
	if(result)
	  message(FATAL_ERROR "Build step for googletest failed: ${result}")
	endif()
	# Prevent overriding the parent project's compiler/linker
	# settings on Windows
	set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

	# Add googletest directly to our build. This defines
	# the gtest and gtest_main targets.
	add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
					 ${CMAKE_BINARY_DIR}/googletest-build)
endif()

#group gtest targets in a folder to avoid crowding the IDE
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_target_properties(gtest PROPERTIES FOLDER thirdparty)
set_target_properties(gtest_main PROPERTIES FOLDER thirdparty)
set_target_properties(gmock PROPERTIES FOLDER thirdparty)
set_target_properties(gmock_main PROPERTIES FOLDER thirdparty)