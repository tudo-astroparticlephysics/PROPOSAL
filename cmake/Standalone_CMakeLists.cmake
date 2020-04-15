CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

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

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

### end full RPATH

SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

SET(DEBUG OFF)
OPTION (DEBUG "DEBUG" OFF)
IF(DEBUG)
        SET( CMAKE_CXX_FLAGS "-g -O0" )
ENDIF()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Some additional options
OPTION(ADD_PYTHON "Choose to compile the python wrapper library" ON)
OPTION(ADD_ROOT "Choose to compile ROOT examples." OFF)
OPTION(ADD_PERFORMANCE_TEST "Choose to compile the performace test source." OFF)
OPTION(ADD_TESTS "Build all unittests." OFF)
OPTION(ADD_CPPEXAMPLE "Choose to compile Cpp example." ON)

#################################################################
#################           python      #########################
#################################################################

IF(ADD_PYTHON)
    MESSAGE(STATUS "Enabled to build the python wrapper library.")
    FIND_PACKAGE( PythonLibs REQUIRED )
    INCLUDE_DIRECTORIES( ${PYTHON_INCLUDE_DIRS} )
	add_subdirectory("vendor/pybind/pybind11")
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

add_subdirectory("vendor/log4cplus/log4cplus")
ADD_DEFINITIONS(-DLOG4CPLUS_SUPPORT=1)
ADD_DEFINITIONS(-DLOG4CPLUS_CONFIG=\"${PROJECT_SOURCE_DIR}/resources/log4cplus.conf\")


#################################################################
#################           Libraries    ########################
#################################################################

SET(LIBRARYS_TO_LINK ${LIBRARYS_TO_LINK} ${CMAKE_THREAD_LIBS_INIT})
SET(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/public/")

INCLUDE_DIRECTORIES(${PROJECT_INCLUDE_DIR})
MESSAGE(STATUS ${PROJECT_INCLUDE_DIR})

FILE(GLOB_RECURSE SRC_FILES ${PROJECT_SOURCE_DIR}/private/PROPOSAL/*)
ADD_LIBRARY(PROPOSAL SHARED ${SRC_FILES})

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    SET_TARGET_PROPERTIES(PROPOSAL PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -Wall -Wextra -Wnarrowing -Wpedantic -fdiagnostics-show-option -Wno-format-security")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    SET_TARGET_PROPERTIES(PROPOSAL PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -Wall -Wextra -Wnarrowing -Wpedantic -fdiagnostics-show-option")
endif()
# The following warnings are silenced, because they arise in the dependencies:
# -Wno-c++11-long-long : "long long" only occurs in ROOT 5

TARGET_LINK_LIBRARIES(PROPOSAL log4cplus ${LIBRARYS_TO_LINK})

IF(IS_SYMLINK ${CMAKE_BINARY_DIR}/resources)
    # Do nothing
ELSE()
    EXECUTE_PROCESS(COMMAND ln -sv ${CMAKE_SOURCE_DIR}/resources ${CMAKE_BINARY_DIR}/resources OUTPUT_VARIABLE link_msg OUTPUT_STRIP_TRAILING_WHITESPACE)
    MESSAGE(STATUS "Symlink to resources created:")
    MESSAGE(STATUS "  ${link_msg}")
ENDIF()

INSTALL(TARGETS PROPOSAL DESTINATION lib)

FILE(GLOB_RECURSE INC_FILES ${PROJECT_SOURCE_DIR}/public/PROPOSAL/*)
# INSTALL(FILES ${INC_FILES} DESTINATION include/PROPOSAL)

foreach(INC_FILE ${INC_FILES})
    file(RELATIVE_PATH REL_FILE ${PROJECT_SOURCE_DIR}/public ${INC_FILE})
    GET_FILENAME_COMPONENT(DIR ${REL_FILE} DIRECTORY)
    INSTALL(FILES "public/${REL_FILE}" DESTINATION include/${DIR})
endforeach()

#################################################################
#################           Executables        ##################
#################################################################

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

IF(ADD_TESTS)
    MESSAGE(STATUS "Building tests enabled.")
	add_subdirectory("vendor/google/googletest" "extern/googletest" EXCLUDE_FROM_ALL)
	mark_as_advanced(
		BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
		gmock_build_tests gtest_build_samples gtest_build_tests
		gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
	)


	set_target_properties(gtest PROPERTIES FOLDER extern)
	set_target_properties(gtest_main PROPERTIES FOLDER extern)
	set_target_properties(gmock PROPERTIES FOLDER extern)
	set_target_properties(gmock_main PROPERTIES FOLDER extern)

	macro(package_add_test TESTNAME)
		add_executable(${TESTNAME} ${ARGN})
		target_link_libraries(${TESTNAME} gtest gmock gtest_main PROPOSAL)
		add_test(NAME ${TESTNAME} COMMAND ${TESTNAME})
		set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
	endmacro()


    ENABLE_TESTING()

    EXECUTE_PROCESS(COMMAND mkdir -p ${PROPOSAL_BINARY_DIR}/bin/ OUTPUT_VARIABLE _output OUTPUT_STRIP_TRAILING_WHITESPACE)

    #create tar directory with "tar -czvf TestFiles.tar.Z TestFiles/" and put it in Test directory
    EXECUTE_PROCESS(COMMAND  tar -xvf ${PROJECT_SOURCE_DIR}/tests/TestFiles.tar.gz -C ${PROPOSAL_BINARY_DIR}/bin/
                    OUTPUT_VARIABLE _output OUTPUT_STRIP_TRAILING_WHITESPACE)

    package_add_test(UnitTest_Utility tests/Utility_TEST.cxx)
    package_add_test(UnitTest_Scattering tests/Scattering_TEST.cxx)
    package_add_test(UnitTest_Photonuclear tests/Photonuclear_TEST.cxx)
    package_add_test(UnitTest_Integral tests/Integral_TEST.cxx)
    package_add_test(UnitTest_Interpolant tests/Interpolant_TEST.cxx)
    package_add_test(UnitTest_Bremsstrahlung tests/Bremsstrahlung_TEST.cxx)
    package_add_test(UnitTest_Compton tests/Compton_TEST.cxx)
    package_add_test(UnitTest_Epairproduction tests/Epairproduction_TEST.cxx)
    package_add_test(UnitTest_Mupairproduction tests/Mupairproduction_TEST.cxx)
    package_add_test(UnitTest_WeakInteraction tests/WeakInteraction_TEST.cxx)
    package_add_test(UnitTest_Annihilation tests/Annihilation_TEST.cxx)
    package_add_test(UnitTest_Ionization tests/Ionization_TEST.cxx)
    package_add_test(UnitTest_PhotoPair tests/PhotoPair_TEST.cxx)
    package_add_test(UnitTest_Medium tests/Medium_TEST.cxx)
    package_add_test(UnitTest_Particle tests/Particle_TEST.cxx)
    package_add_test(UnitTest_ParticleDef tests/ParticleDef_TEST.cxx)
    package_add_test(UnitTest_DecayChannel tests/DecayChannel_TEST.cxx)
    package_add_test(UnitTest_DecayTable tests/DecayTable_TEST.cxx)
    package_add_test(UnitTest_EnergyCutSettings tests/EnergyCutSettings_TEST.cxx)
    package_add_test(UnitTest_ContinuousRandomization tests/ContinuousRandomization_TEST.cxx)
    package_add_test(UnitTest_Geometry tests/Geometry_TEST.cxx)
    package_add_test(UnitTest_Vector3D tests/Vector3D_TEST.cxx)
    package_add_test(UnitTest_Propagation tests/Propagation_TEST.cxx)
    package_add_test(UnitTest_Sector tests/Sector_TEST.cxx)
    package_add_test(UnitTest_MathMethods tests/MathMethods_TEST.cxx)
    package_add_test(UnitTest_Spline tests/Spline_TEST.cxx)
    package_add_test(UnitTest_Density tests/Density_distribution_TEST.cxx)
    package_add_test(UnitTest_Displacement tests/Displacement_TEST.cxx)
    package_add_test(UnitTest_Time tests/Time_TEST.cxx)
    package_add_test(UnitTest_Interaction tests/Interaction_TEST.cxx)
    package_add_test(UnitTest_UtilityIntegral tests/UtilityIntegral_TEST.cxx)
    package_add_test(UnitTest_UtilityInterpolant tests/UtilityInterpolant_TEST.cxx)

else (ADD_TESTS)
    MESSAGE(STATUS "No tests will be build.")
ENDIF()

ADD_SUBDIRECTORY( doc )


IF(ADD_PYTHON)
	FILE(GLOB_RECURSE PYTHON_SRC_FILES
		${PROJECT_SOURCE_DIR}/private/Interface/python/*)
	INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/public/Interface/python )
	PYBIND11_ADD_MODULE(pyPROPOSAL SHARED ${PYTHON_SRC_FILES})
	TARGET_LINK_LIBRARIES(pyPROPOSAL PRIVATE PROPOSAL)
	SET_TARGET_PROPERTIES(pyPROPOSAL PROPERTIES PREFIX "" SUFFIX ".so"
		COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")
	INSTALL(TARGETS pyPROPOSAL DESTINATION lib)
ENDIF(ADD_PYTHON)

# uninstall target
CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
