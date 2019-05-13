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

#################################################################
#################           python      #########################
#################################################################

IF(ADD_PYTHON)
	MESSAGE(STATUS "Enabled to build the python wrapper library.")
	FIND_PACKAGE( PythonLibs 2.7 REQUIRED )

	IF(PYTHONLIBS_FOUND)
		# TODO(mario): Find a better way to search for pybind11!
		# Preferably with a FindPybind11.make.
		# I haven't found a reliable method yet.

		INCLUDE_DIRECTORIES( ${PYTHON_INCLUDE_DIRS} )

		FIND_PACKAGE(pybind11 CONFIG)

		IF(pybind11_FOUND)
			MESSAGE(STATUS "Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIR}")
			INCLUDE_DIRECTORIES( ${pybind11_INCLUDE_DIR} )
		ELSE(pybind11_FOUND)
			MESSAGE(STATUS "Clone pybind11")

			execute_process(
				COMMAND git clone https://github.com/pybind/pybind11 --branch v2.2.4 ${CMAKE_CURRENT_SOURCE_DIR}/pybind11
				OUTPUT_VARIABLE     git_output
				ERROR_VARIABLE		git_error
			)

			add_subdirectory(pybind11)

			SET(pybind11_FOUND TRUE)
			MESSAGE(STATUS "Clone pybind11: ${git_error}")
		ENDIF(pybind11_FOUND)

	ENDIF(PYTHONLIBS_FOUND)
ELSE(ADD_PYTHON)
	MESSAGE(STATUS "No python wrapper library will be build.")
ENDIF(ADD_PYTHON)

#################################################################
#################           boost       #########################
#################################################################

FIND_PACKAGE( Boost REQUIRED )
SET(LIBRARYS_TO_LINK ${LIBRARYS_TO_LINK} ${Boost_LIBRARIES})

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
		# MESSAGE(STATUS "ROOT examples will not be builded.")
		MESSAGE(STATUS "No ROOT Output is available.")
		MESSAGE(STATUS "Make sure you have ROOT installed and ROOTSYS is set.")
	ENDIF(ROOT_FOUND)
ENDIF(ADD_ROOT)


#################################################################
#################           log4cplus     #######################
#################################################################

INCLUDE(cmake/FindLog4cplus.cmake)

IF(LOG4CPLUS_FOUND)
    ADD_DEFINITIONS(-DLOG4CPLUS_SUPPORT=1)
	ADD_DEFINITIONS(-DLOG4CPLUS_CONFIG=\"${PROJECT_SOURCE_DIR}/resources/log4cplus.conf\")
    SET(LIBRARYS_TO_LINK ${LIBRARYS_TO_LINK} ${LOG4CPLUS_LIBRARIES})
ELSE(LOG4CPLUS_FOUND)
    ADD_DEFINITIONS(-DLOG4CPLUS_SUPPORT=0)
    MESSAGE(STATUS "No logging will be done.")
ENDIF(LOG4CPLUS_FOUND)


#################################################################
#################           GTest       #########################
#################################################################

if (ADD_TESTS)
    MESSAGE(STATUS "Building tests enabled.")
	INCLUDE(cmake/gtest.cmake)
	SET(gtest_LIBRARIES libgtest)
 	ENABLE_TESTING()
else (ADD_TESTS)
    MESSAGE(STATUS "No tests will be build.")
endif (ADD_TESTS)

#################################################################
#################           Libraries    ########################
#################################################################

SET(LIBRARYS_TO_LINK ${LIBRARYS_TO_LINK} ${CMAKE_THREAD_LIBS_INIT})

INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/public" "${PROJECT_SOURCE_DIR}" ${LOG4CPLUS_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

# FILE(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/private/PROPOSAL/*.cxx)
FILE(GLOB_RECURSE SRC_FILES ${PROJECT_SOURCE_DIR}/private/PROPOSAL/*)
ADD_LIBRARY(PROPOSAL SHARED ${SRC_FILES})

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	SET_TARGET_PROPERTIES(PROPOSAL PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -Wall -Wextra -Wnarrowing -Wpedantic -fdiagnostics-show-option -Wno-format-security")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	SET_TARGET_PROPERTIES(PROPOSAL PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -Wall -Wextra -Wnarrowing -Wpedantic -fdiagnostics-show-option")
endif()
# The following warnings are silenced, because they arise in the dependencies:
# -Wno-c++11-long-long : "long long" only occurs in ROOT 5
# -Wno-variadic-macros : only in output.h; has something todo with BOOST
# -Wno-gnu-zero-variadic-macro-arguments : this is also a warning comming from BOOST occuring e.g. in log_fatal()

TARGET_LINK_LIBRARIES(PROPOSAL ${LIBRARYS_TO_LINK})

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

#################################################################
#################           Executables        ##################
#################################################################

ADD_EXECUTABLE(WriteSectorsFromDomList
        private/test/WriteSectorsFromDomList.cxx
)
TARGET_LINK_LIBRARIES(WriteSectorsFromDomList PROPOSAL)

ADD_EXECUTABLE(example
        private/test/example.cxx
)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	SET_TARGET_PROPERTIES(example PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -Wall -Wextra -Wnarrowing -Wpedantic -fdiagnostics-show-option")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	SET_TARGET_PROPERTIES(example PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -Wall -Wextra -Wnarrowing -Wpedantic -fdiagnostics-show-option")
endif()

TARGET_LINK_LIBRARIES(example PROPOSAL)

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
  EXECUTE_PROCESS(COMMAND mkdir -p ${PROPOSAL_BINARY_DIR}/bin/ OUTPUT_VARIABLE _output OUTPUT_STRIP_TRAILING_WHITESPACE)

  #create tar directory with "tar -czvf TestFiles.tar.Z TestFiles/" and put it in Test directory
  EXECUTE_PROCESS(COMMAND  tar -xvf ${PROJECT_SOURCE_DIR}/tests/TestFiles.tar.gz -C ${PROPOSAL_BINARY_DIR}/bin/
                    OUTPUT_VARIABLE _output OUTPUT_STRIP_TRAILING_WHITESPACE)

  ADD_EXECUTABLE(UnitTest_Utility tests/Utility_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_Scattering tests/Scattering_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_Photonuclear tests/Photonuclear_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_Integral tests/Integral_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_Interpolant tests/Interpolant_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_Bremsstrahlung tests/Bremsstrahlung_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_Epairproduction tests/Epairproduction_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_Mupairproduction tests/Mupairproduction_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_WeakInteraction tests/WeakInteraction_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_Ionization tests/Ionization_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_Medium tests/Medium_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_Particle tests/Particle_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_ParticleDef tests/ParticleDef_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_DecayChannel tests/DecayChannel_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_DecayTable tests/DecayTable_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_EnergyCutSettings tests/EnergyCutSettings_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_ContinuousRandomization tests/ContinuousRandomization_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_Geometry tests/Geometry_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_Vector3D tests/Vector3D_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_Propagation tests/Propagation_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_Sector tests/Sector_TEST.cxx)

  TARGET_LINK_LIBRARIES(UnitTest_Utility PROPOSAL ${gtest_LIBRARIES})
  TARGET_LINK_LIBRARIES(UnitTest_Scattering PROPOSAL ${gtest_LIBRARIES})
  TARGET_LINK_LIBRARIES(UnitTest_Integral PROPOSAL ${gtest_LIBRARIES})
  TARGET_LINK_LIBRARIES(UnitTest_Interpolant PROPOSAL ${gtest_LIBRARIES})
  TARGET_LINK_LIBRARIES(UnitTest_Ionization PROPOSAL ${gtest_LIBRARIES})
  TARGET_LINK_LIBRARIES(UnitTest_Bremsstrahlung PROPOSAL ${gtest_LIBRARIES})
  TARGET_LINK_LIBRARIES(UnitTest_Epairproduction PROPOSAL ${gtest_LIBRARIES})
  TARGET_LINK_LIBRARIES(UnitTest_Mupairproduction PROPOSAL ${gtest_LIBRARIES})
  TARGET_LINK_LIBRARIES(UnitTest_WeakInteraction PROPOSAL ${gtest_LIBRARIES})
  TARGET_LINK_LIBRARIES(UnitTest_Photonuclear PROPOSAL ${gtest_LIBRARIES})
  TARGET_LINK_LIBRARIES(UnitTest_Medium PROPOSAL ${gtest_LIBRARIES})
  TARGET_LINK_LIBRARIES(UnitTest_Particle PROPOSAL ${gtest_LIBRARIES})
  TARGET_LINK_LIBRARIES(UnitTest_ParticleDef PROPOSAL ${gtest_LIBRARIES})
  TARGET_LINK_LIBRARIES(UnitTest_DecayChannel PROPOSAL ${gtest_LIBRARIES})
  TARGET_LINK_LIBRARIES(UnitTest_DecayTable PROPOSAL ${gtest_LIBRARIES})
  TARGET_LINK_LIBRARIES(UnitTest_EnergyCutSettings PROPOSAL ${gtest_LIBRARIES})
  TARGET_LINK_LIBRARIES(UnitTest_ContinuousRandomization PROPOSAL ${gtest_LIBRARIES})
  TARGET_LINK_LIBRARIES(UnitTest_Geometry PROPOSAL ${gtest_LIBRARIES})
  TARGET_LINK_LIBRARIES(UnitTest_Vector3D PROPOSAL ${gtest_LIBRARIES})
  TARGET_LINK_LIBRARIES(UnitTest_Propagation PROPOSAL ${gtest_LIBRARIES})
  TARGET_LINK_LIBRARIES(UnitTest_Sector PROPOSAL ${gtest_LIBRARIES})

  ADD_TEST(UnitTest_Utility bin/UnitTest_Utility)
  ADD_TEST(UnitTest_Scattering bin/UnitTest_Scattering)
  ADD_TEST(UnitTest_ContinuousRandomization bin/UnitTest_ContinuousRandomization)
  ADD_TEST(UnitTest_Integral bin/UnitTest_Integral)
  ADD_TEST(UnitTest_Medium bin/UnitTest_Medium)
  ADD_TEST(UnitTest_Particle bin/UnitTest_Particle)
  ADD_TEST(UnitTest_ParticleDef bin/UnitTest_ParticleDef)
  ADD_TEST(UnitTest_DecayChannel bin/UnitTest_DecayChannel)
  ADD_TEST(UnitTest_DecayTable bin/UnitTest_DecayTable)
  ADD_TEST(UnitTest_EnergyCutSettings bin/UnitTest_EnergyCutSettings)
  ADD_TEST(UnitTest_Interpolant bin/UnitTest_Interpolant)
  ADD_TEST(UnitTest_Epairproduction bin/UnitTest_Epairproduction)
  ADD_TEST(UnitTest_Mupairproduction bin/UnitTest_Mupairproduction)
  ADD_TEST(UnitTest_WeakInteraction bin/UnitTest_WeakInteraction)
  ADD_TEST(UnitTest_Ionization bin/UnitTest_Ionization)
  ADD_TEST(UnitTest_Bremsstrahlung bin/UnitTest_Bremsstrahlung)
  ADD_TEST(UnitTest_Photonuclear bin/UnitTest_Photonuclear)
  ADD_TEST(UnitTest_Geometry bin/UnitTest_Geometry)
  ADD_TEST(UnitTest_Vector3D bin/UnitTest_Vector3D)
  ADD_TEST(UnitTest_Propagation bin/UnitTest_Propagation)
  ADD_TEST(UnitTest_Sector bin/UnitTest_Sector)

ENDIF()

ADD_SUBDIRECTORY( doc )


IF(ADD_PYTHON)
	IF(PYTHONLIBS_FOUND AND pybind11_FOUND)
		PYBIND11_ADD_MODULE(pyPROPOSAL SHARED private/python/pybindings.cxx)
		TARGET_LINK_LIBRARIES(pyPROPOSAL PRIVATE PROPOSAL)
		SET_TARGET_PROPERTIES(pyPROPOSAL PROPERTIES PREFIX "" SUFFIX ".so" COMPILE_FLAGS "${CMAKE_CXX_FLAGS}")
		INSTALL(TARGETS pyPROPOSAL DESTINATION lib)
	ENDIF(PYTHONLIBS_FOUND AND pybind11_FOUND)
ENDIF(ADD_PYTHON)

# uninstall target
CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

