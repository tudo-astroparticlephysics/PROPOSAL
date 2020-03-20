cmake_minimum_required(VERSION 3.8)

IF(APPLE)
    # In newer version of cmake this will be the default
    SET(CMAKE_MACOSX_RPATH 1)
ENDIF(APPLE)

# sets standard installtion paths
include(GNUInstallDirs)

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
OPTION(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
list(APPEND CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${CMAKE_INSTALL_PREFIX}/lib64")
message(STATUS "${CMAKE_INSTALL_RPATH}")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
OPTION(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
### end full RPATH

# set a default build type (Debug for git checkout else Release)
include(cmake/BuildType.cmake)

# Some additional options
OPTION(ADD_PYTHON "Choose to compile the python wrapper library" ON)
OPTION(ADD_ROOT "Choose to compile ROOT examples." OFF)
OPTION(ADD_PERFORMANCE_TEST "Choose to compile the performace test source." OFF)
OPTION(ADD_CPPEXAMPLE "Choose to compile Cpp example." ON)


#################################################################
####################       PROPOSAL      ########################
#################################################################
file(GLOB_RECURSE SRC_FILES ${PROJECT_SOURCE_DIR}/private/PROPOSAL/*)
add_library(PROPOSAL SHARED ${SRC_FILES})
add_library(PROPOSAL::PROPOSAL ALIAS PROPOSAL)
target_compile_features(PROPOSAL PUBLIC cxx_std_11)
set_target_properties(PROPOSAL PROPERTIES CXX_EXTENSIONS OFF)
target_include_directories(
    PROPOSAL PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/public>
    $<INSTALL_INTERFACE:include>
)
target_compile_options(PROPOSAL PRIVATE -Wall -Wextra -Wnarrowing -Wpedantic -fdiagnostics-show-option -Wno-format-security)
install(
    TARGETS PROPOSAL
    EXPORT PROPOSALTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# input version from the project call, so the library knows its own version
configure_file(
    "${PROJECT_SOURCE_DIR}/public/PROPOSAL/version.h.in"
    "${PROJECT_BINARY_DIR}/include/PROPOSAL/version.h"
)
install(
    FILES ${PROJECT_BINARY_DIR}/include/PROPOSAL/version.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/PROPOSAL
)
target_include_directories(
    PROPOSAL PUBLIC
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# install header files
file(GLOB_RECURSE INC_FILES ${PROJECT_SOURCE_DIR}/public/PROPOSAL/*)
foreach(INC_FILE ${INC_FILES})
    file(RELATIVE_PATH REL_FILE ${PROJECT_SOURCE_DIR}/public ${INC_FILE})
    get_filename_component(DIR ${REL_FILE} DIRECTORY)
    install(FILES "public/${REL_FILE}" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${DIR})
endforeach()


if(NOT IS_SYMLINK ${CMAKE_BINARY_DIR}/resources)
    execute_process(COMMAND ln -sv ${CMAKE_SOURCE_DIR}/resources ${CMAKE_BINARY_DIR}/resources OUTPUT_VARIABLE link_msg OUTPUT_STRIP_TRAILING_WHITESPACE)
    message(STATUS "Symlink to resources created:")
    message(STATUS "  ${link_msg}")
endif()


#################################################################
#################           log4cplus     #######################
#################################################################

option(LOG4CPLUS_BUILD_TESTING "Build unit tests for log4cplus" OFF)
add_subdirectory("vendor/log4cplus/log4cplus")
target_link_libraries(PROPOSAL PRIVATE log4cplus)
target_compile_definitions(PROPOSAL PRIVATE -DLOG4CPLUS_SUPPORT=1)

#################################################################
#################           Executables        ##################
#################################################################

IF(ADD_CPPEXAMPLE)
    add_executable(example private/test/example.cxx)
    target_compile_options(example PRIVATE -Wall -Wextra -Wnarrowing -Wpedantic -fdiagnostics-show-option)
    target_link_libraries(example PRIVATE PROPOSAL)
ENDIF(ADD_CPPEXAMPLE)


IF(ADD_PERFORMANCE_TEST)
    add_executable(performance_test private/test/performance_test.cxx)
    target_compile_options(performance_test PRIVATE -Wall -Wextra -Wnarrowing -Wpedantic -fdiagnostics-show-option)
    target_link_libraries(performance_test PRIVATE PROPOSAL)
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

#################################################################
#################           Tests        ########################
#################################################################
add_subdirectory(doc)


#################################################################
#################           python      #########################
#################################################################
IF(ADD_PYTHON)
    message(STATUS "Building the python wrapper library.")
    find_package(PythonLibs REQUIRED)
    add_subdirectory("vendor/pybind/pybind11")
    file(GLOB_RECURSE PYTHON_SRC_FILES ${PROJECT_SOURCE_DIR}/private/Interface/python/*)

    pybind11_add_module(pyproposal SHARED ${PYTHON_SRC_FILES})
    set_target_properties(pyproposal PROPERTIES OUTPUT_NAME proposal)
    target_include_directories(pyproposal PRIVATE ${PYTHON_INCLUDE_DIRS})
    target_include_directories(pyproposal PRIVATE ${PROJECT_SOURCE_DIR}/public/Interface/python )
    target_link_libraries(pyproposal PRIVATE PROPOSAL)
    target_compile_options(pyproposal PRIVATE -fvisibility=hidden)
    install(TARGETS pyproposal EXPORT PROPOSALTargets DESTINATION ${CMAKE_INSTALL_LIBDIR})
ELSE(ADD_PYTHON)
    MESSAGE(STATUS "No python wrapper library will be build.")
ENDIF(ADD_PYTHON)


#################################################################
#################           ROOT        #########################
#################################################################

IF(ADD_ROOT)
    message(STATUS "Enabling ROOT support.")

    # Load some basic macros which are needed later on
    list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
    find_package(ROOT REQUIRED)

    target_compile_definitions(PROPOSAL PRIVATE -DROOT_SUPPORT=1)
    target_include_directories(PROPOSAL PRIVATE ${ROOT_INCLUDE_DIR})
    target_link_libraries(PROPOSAL PRIVATE ${ROOT_LIBRARIES})
ENDIF(ADD_ROOT)

# uninstall target
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY
)
add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
)



#################################################################
#################         INSTALLATION        ###################
#################################################################
install(
    EXPORT PROPOSALTargets
    FILE PROPOSALConfig.cmake
    NAMESPACE PROPOSAL::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PROPOSAL
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "PROPOSALConfigVersion.cmake"
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY SameMajorVersion
)
install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/PROPOSALConfigVersion.cmake"
    DESTINATION lib/cmake/PROPOSAL
)
