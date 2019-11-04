CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF(APPLE)
    # In newer version of cmake this will be the default
    SET(CMAKE_MACOSX_RPATH 1)
ENDIF(APPLE)

### full RPATH
### copied from https://cmake.org/Wiki/CMake_RPATH_handling
### set the RPATH so that for using PROPOSAL in python
### the (DY)LD_LIBRARY_PATH don't have to be be set 
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

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")
### end full RPATH

# Some additional options
OPTION(ADD_PYTHON "Choose to compile the python wrapper library" ON)
OPTION(ADD_ROOT "Choose to compile ROOT examples." OFF)
OPTION(ADD_PERFORMANCE_TEST "Choose to compile the performace test source." OFF)
OPTION(ADD_CPPEXAMPLE "Choose to compile Cpp example." ON)

#################################################################
#################           python      #########################
#################################################################

IF(ADD_PYTHON)
    MESSAGE(STATUS "Enabled to build the python wrapper library.")
    FIND_PACKAGE(PythonLibs 2.7 REQUIRED)

    IF(PYTHONLIBS_FOUND)
		add_subdirectory(extern/pybind11)
    ENDIF(PYTHONLIBS_FOUND)
ELSE(ADD_PYTHON)
    MESSAGE(STATUS "No python wrapper library will be build.")
ENDIF(ADD_PYTHON)

#################################################################
#################           ROOT        #########################
#################################################################

IF(ADD_ROOT)
    MESSAGE(STATUS "Enabled ROOT support.")
    # Load some basic macros which are needed later on
    INCLUDE(cmake/FindROOT_new.cmake)

    #if ROOT is found ROOT files with ROOT trees can be written
    IF(ROOT_FOUND)
        ADD_DEFINITIONS(-DROOT_SUPPORT=1)

        INCLUDE_DIRECTORIES(${ROOT_INCLUDE_DIR})
        SET(LIBRARYS_TO_LINK ${LIBRARYS_TO_LINK} ${ROOT_LIBRARIES})

    ELSE(ROOT_FOUND)
        ADD_DEFINITIONS(-DROOT_SUPPORT=0)

        MESSAGE(STATUS "ROOT not found...")
        MESSAGE(STATUS "No ROOT Output is available.")
        MESSAGE(STATUS "Make sure you have ROOT installed and ROOTSYS is set.")
    ENDIF(ROOT_FOUND)
ENDIF(ADD_ROOT)


#################################################################
#################           log4cplus     #######################
#################################################################

add_subdirectory(extern/log4cplus)


#################################################################
#################           Libraries    ########################
#################################################################

FILE(GLOB_RECURSE SRC_FILES ${PROJECT_SOURCE_DIR}/private/PROPOSAL/*)
add_library(PROPOSAL SHARED ${SRC_FILES})
target_include_directories(PROPOSAL PUBLIC ${PROJECT_SOURCE_DIR}/public)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    SET_TARGET_PROPERTIES(PROPOSAL PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -Wall -Wextra -Wnarrowing -Wpedantic -fdiagnostics-show-option -Wno-format-security")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    SET_TARGET_PROPERTIES(PROPOSAL PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -Wall -Wextra -Wnarrowing -Wpedantic -fdiagnostics-show-option")
endif()
# The following warnings are silenced, because they arise in the dependencies:
# -Wno-c++11-long-long : "long long" only occurs in ROOT 5

TARGET_LINK_LIBRARIES(PROPOSAL log4cplus)

IF(IS_SYMLINK ${CMAKE_BINARY_DIR}/resources)
    # Do nothing
ELSE()
    EXECUTE_PROCESS(COMMAND ln -sv ${CMAKE_SOURCE_DIR}/resources ${CMAKE_BINARY_DIR}/resources OUTPUT_VARIABLE link_msg OUTPUT_STRIP_TRAILING_WHITESPACE)
    MESSAGE(STATUS "Symlink to resources created:")
    MESSAGE(STATUS "  ${link_msg}")
ENDIF()

INSTALL(TARGETS PROPOSAL DESTINATION lib)

FILE(GLOB_RECURSE INC_FILES ${PROJECT_SOURCE_DIR}/public/PROPOSAL/*.h)
# INSTALL(FILES ${INC_FILES} DESTINATION include/PROPOSAL)

foreach(INC_FILE ${INC_FILES})
    file(RELATIVE_PATH REL_FILE ${PROJECT_SOURCE_DIR}/public ${INC_FILE})
    GET_FILENAME_COMPONENT(DIR "${REL_FILE}" DIRECTORY)
    INSTALL(FILES "public/${REL_FILE}" DESTINATION include/${DIR})
endforeach()

ADD_EXECUTABLE(WriteSectorsFromDomList
    private/test/WriteSectorsFromDomList.cxx
)
TARGET_LINK_LIBRARIES(WriteSectorsFromDomList PROPOSAL)

IF(ADD_CPPEXAMPLE)
    ADD_EXECUTABLE(example
            private/test/example.cxx
    )
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        SET_TARGET_PROPERTIES(example PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -Wall -Wextra -Wnarrowing -Wpedantic -fdiagnostics-show-option")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        SET_TARGET_PROPERTIES(example PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -Wall -Wextra -Wnarrowing -Wpedantic -fdiagnostics-show-option")
endif()

TARGET_LINK_LIBRARIES(example PROPOSAL)
ENDIF(ADD_CPPEXAMPLE)


IF(ADD_PERFORMANCE_TEST)
    ADD_EXECUTABLE(performance_test
        private/test/performance_test.cxx
    )
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        SET_TARGET_PROPERTIES(performance_test PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -Wall -Wextra -Wnarrowing -Wpedantic -fdiagnostics-show-option")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        SET_TARGET_PROPERTIES(performance_test PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -Wall -Wextra -Wnarrowing -Wpedantic -fdiagnostics-show-option")
    endif()
    TARGET_LINK_LIBRARIES(performance_test PROPOSAL ${LIBRARYS_TO_LINK_PERFORMANCE_TEST})
ENDIF(ADD_PERFORMANCE_TEST)

#################################################################
#################           Tests        ########################
#################################################################

# only include tests if PROPOSAL is not build as a dependency
if(CMAKE_PROJECT_NAME STREQUAL PROPOSAL)
	INCLUDE(CTest)
endif()
if(CMAKE_PROJECT_NAME STREQUAL PROPOSAL AND BUILD_TESTING)
	add_subdirectory(tests)
endif()

ADD_SUBDIRECTORY( doc )


IF(ADD_PYTHON)
    IF(PYTHONLIBS_FOUND)
        FILE(GLOB_RECURSE PYTHON_SRC_FILES ${PROJECT_SOURCE_DIR}/private/python/*)
        PYBIND11_ADD_MODULE(pyPROPOSAL SHARED ${PYTHON_SRC_FILES})
		TARGET_INCLUDE_DIRECTORIES(pyPROPOSAL PUBLIC ${PROJECT_SOURCE_DIR}/public/python)
		TARGET_LINK_LIBRARIES(pyPROPOSAL PUBLIC PROPOSAL)
        # SET_TARGET_PROPERTIES(pyPROPOSAL PROPERTIES PREFIX "" SUFFIX ".so"
        #     COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")
        INSTALL(TARGETS pyPROPOSAL DESTINATION lib)
    ENDIF(PYTHONLIBS_FOUND)
ENDIF(ADD_PYTHON)

# uninstall target
CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
