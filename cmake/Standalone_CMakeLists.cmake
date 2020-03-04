CMAKE_MINIMUM_REQUIRED(VERSION 3.8)

IF(APPLE)
    # In newer version of cmake this will be the default
    SET(CMAKE_MACOSX_RPATH 1)
ENDIF(APPLE)

### full RPATH
### copied from https://cmake.org/Wiki/CMake_RPATH_handling
### set the RPATH so that for using PROPOSAL in python
### the DYLD_LIBRARY_PATH must not be set in the bashrc
### But for using PROPOSAL as c-Library, this path still
### has to be set

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
### end full RPATH

SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# set a default build type (Debug for git checkout else Release)
include(cmake/BuildType.cmake)

# Some additional options
OPTION(ADD_PYTHON "Choose to compile the python wrapper library" ON)
OPTION(ADD_ROOT "Choose to compile ROOT examples." OFF)
OPTION(ADD_PERFORMANCE_TEST "Choose to compile the performace test source." OFF)
OPTION(ADD_CPPEXAMPLE "Choose to compile Cpp example." ON)

#################################################################
#################           ROOT        #########################
#################################################################

IF(ADD_ROOT)
    MESSAGE(STATUS "Enabling ROOT support.")
    # Load some basic macros which are needed later on
    INCLUDE(cmake/FindROOT_new.cmake)

    #if ROOT is found ROOT files with ROOT trees can be written
    IF(ROOT_FOUND)
        ADD_DEFINITIONS(-DROOT_SUPPORT=1)

        INCLUDE_DIRECTORIES(${ROOT_INCLUDE_DIR})
        SET(LIBRARIES_TO_LINK ${LIBRARIES_TO_LINK} ${ROOT_LIBRARIES})

    ELSE(ROOT_FOUND)
        ADD_DEFINITIONS(-DROOT_SUPPORT=0)

		MESSAGE(SEND_ERROR  "ROOT not found...")
		MESSAGE(FATAL_ERROR "Make sure you have ROOT installed and ROOTSYS is set.")
    ENDIF(ROOT_FOUND)
ENDIF(ADD_ROOT)


#################################################################
#################           log4cplus     #######################
#################################################################

option(LOG4CPLUS_BUILD_TESTING "Build unit tests for log4cplus" OFF)
add_subdirectory("vendor/log4cplus/log4cplus")
ADD_DEFINITIONS(-DLOG4CPLUS_SUPPORT=1)
ADD_DEFINITIONS(-DLOG4CPLUS_CONFIG=\"${PROJECT_SOURCE_DIR}/resources/log4cplus.conf\")


#################################################################
#################           Libraries    ########################
#################################################################

SET(LIBRARIES_TO_LINK ${LIBRARIES_TO_LINK} ${CMAKE_THREAD_LIBS_INIT})
SET(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/public/")


FILE(GLOB_RECURSE SRC_FILES ${PROJECT_SOURCE_DIR}/private/PROPOSAL/*)

add_library(PROPOSAL SHARED ${SRC_FILES})
target_compile_features(PROPOSAL PUBLIC cxx_std_11)
set_target_properties(PROPOSAL PROPERTIES CXX_EXTENSIONS OFF)
target_include_directories(PROPOSAL PUBLIC ${PROJECT_INCLUDE_DIR})
target_compile_options(PROPOSAL PRIVATE -Wall -Wextra -Wnarrowing -Wpedantic -fdiagnostics-show-option -Wno-format-security)
target_link_libraries(PROPOSAL log4cplus ${LIBRARIES_TO_LINK})

IF(IS_SYMLINK ${CMAKE_BINARY_DIR}/resources)
    # Do nothing
ELSE()
    execute_process(COMMAND ln -sv ${CMAKE_SOURCE_DIR}/resources ${CMAKE_BINARY_DIR}/resources OUTPUT_VARIABLE link_msg OUTPUT_STRIP_TRAILING_WHITESPACE)
    message(STATUS "Symlink to resources created:")
    message(STATUS "  ${link_msg}")
ENDIF()

install(TARGETS PROPOSAL LIBRARY DESTINATION lib)

FILE(GLOB_RECURSE INC_FILES ${PROJECT_SOURCE_DIR}/public/PROPOSAL/*)

foreach(INC_FILE ${INC_FILES})
    file(RELATIVE_PATH REL_FILE ${PROJECT_SOURCE_DIR}/public ${INC_FILE})
    GET_FILENAME_COMPONENT(DIR ${REL_FILE} DIRECTORY)
    INSTALL(FILES "public/${REL_FILE}" DESTINATION include/${DIR})
endforeach()

#################################################################
#################           Executables        ##################
#################################################################

IF(ADD_CPPEXAMPLE)
    add_executable(example private/test/example.cxx)
	target_compile_options(example PRIVATE -Wall -Wextra -Wnarrowing -Wpedantic -fdiagnostics-show-option)
	target_link_libraries(example PROPOSAL)
ENDIF(ADD_CPPEXAMPLE)


IF(ADD_PERFORMANCE_TEST)
    add_executable(performance_test private/test/performance_test.cxx)
    target_compile_options(performance_test -Wall -Wextra -Wnarrowing -Wpedantic -fdiagnostics-show-option)
    target_link_libraries(performance_test PROPOSAL)
ENDIF(ADD_PERFORMANCE_TEST)

#################################################################
#################           Tests        ########################
#################################################################

if(CMAKE_PROJECT_NAME STREQUAL PROPOSAL)
	include(CTest)
endif()
if(CMAKE_PROJECT_NAME STREQUAL PROPOSAL AND BUILD_TESTING)
    add_subdirectory(tests)
else()
    MESSAGE(STATUS "No tests will be build.")
endif()

ADD_SUBDIRECTORY(doc)


#################################################################
#################           python      #########################
#################################################################
IF(ADD_PYTHON)
    MESSAGE(STATUS "Building the python wrapper library.")
    find_package(PythonLibs REQUIRED)
	add_subdirectory("vendor/pybind/pybind11")

	FILE(GLOB_RECURSE PYTHON_SRC_FILES ${PROJECT_SOURCE_DIR}/private/Interface/python/*)
	PYBIND11_ADD_MODULE(pyPROPOSAL SHARED ${PYTHON_SRC_FILES})
	target_include_directories(pyPROPOSAL PRIVATE ${PYTHON_INCLUDE_DIRS})
	target_include_directories(pyPROPOSAL PRIVATE ${PROJECT_SOURCE_DIR}/public/Interface/python )
	target_link_libraries(pyPROPOSAL PRIVATE PROPOSAL)
	target_compile_options(pyPROPOSAL PRIVATE -fvisibility=hidden)
	INSTALL(TARGETS pyPROPOSAL DESTINATION lib)
ELSE(ADD_PYTHON)
    MESSAGE(STATUS "No python wrapper library will be build.")
ENDIF(ADD_PYTHON)

# uninstall target
CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
