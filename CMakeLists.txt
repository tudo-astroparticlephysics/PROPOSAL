cmake_minimum_required(VERSION 2.6)

#FIND_LIBRARY(LOG4CPLUS_LIB NAMES log4cplus HINTS $LD_LIBRARY_PATH)
#MESSAGE(STATUS "log4cplus lib path: ${LOG4CPLUS_LIB}")
project( PROPOSAL)

# This module defines
# LOG4CPP_INCLUDE_DIR, where to find ptlib.h, etc.
# LOG4CPP_LIBRARIES, the libraries to link against to use pwlib.
# LOG4CPP_FOUND, If false, don't try to use pwlib.

#################################################################
#################           log4cpp     #########################
#################################################################
#FIND_PATH(LOG4CPP_INCLUDE_DIR log4cpp/Category.hh
#  PATHS
#    "$ENV{LOG4CPP}/include"
#    /usr/local/include
#    /usr/include
#)
#
#FIND_LIBRARY(LOG4CPP_LIBRARIES log4cpp
#  PATHS
#    "$ENV{LOG4CPP}/lib"
#    /usr/local/lib
#    /usr/lib
#)
#
#SET(LOG4CPP_FOUND 0)
#IF(LOG4CPP_INCLUDE_DIR)
#  IF(LOG4CPP_LIBRARIES)
#    SET(LOG4CPP_FOUND 1)
#    MESSAGE(STATUS "Found Log4CPP")
#  ENDIF(LOG4CPP_LIBRARIES)
#ENDIF(LOG4CPP_INCLUDE_DIR)

#MARK_AS_ADVANCED(
#  LOG4CPP_INCLUDE_DIR
#  LOG4CPP_LIBRARIES
#)

#################################################################
#################           GTest       #########################
#################################################################
IF(DEFINED ENV{GTEST_PATH})
  enable_testing()
  #1. Path
  #2. Path where the binarys will be copied to
  add_subdirectory($ENV{GTEST_PATH} ${CMAKE_CURRENT_BINARY_DIR}/gtestbin)
  set(DoTesting true)
  MESSAGE(STATUS "Found gtest-path")
ELSE()
    message("GTEST_PATH is not defined.  To build testsuite you have to set GTEST_PATH to gtest root directory")
    message("Remember to add gtest/include to C_INCLUDE_PATH and CPP_INCLUDE_PATH")
ENDIF()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

SET(DEBUG OFF)
option (DEBUG "DEBUG" OFF)
if(DEBUG)
        set( CMAKE_CXX_FLAGS "-g -pg -O0" )
endif()

include_directories("${PROJECT_SOURCE_DIR}/public" ${GTEST_INCLUDE_DIR} ${LOG4CPP_INCLUDE_DIR})

add_library(PROPOSAL
        private/PROPOSAL/Integral.cxx
        private/PROPOSAL/methods.cxx
        private/PROPOSAL/MathModel.cxx
        private/PROPOSAL/Output.cxx
        private/PROPOSAL/Bremsstrahlung.cxx
        private/PROPOSAL/CrossSections.cxx
        private/PROPOSAL/Decay.cxx
        private/PROPOSAL/Epairproduction.cxx
        private/PROPOSAL/Ionization.cxx
        private/PROPOSAL/Photonuclear.cxx
        private/PROPOSAL/Medium.cxx
        private/PROPOSAL/Particle.cxx
        private/PROPOSAL/EnergyCutSettings.cxx
	private/PROPOSAL/Interpolant.cxx
        private/PROPOSAL/StandardNormal.cxx
        private/PROPOSAL/RootFinder.cxx
     #   private/PROPOSAL/ProcessCollection.cxx

	)

set_target_properties(PROPOSAL PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror")# -Wextra -pedantic")

add_executable(PROPOSALtest
        private/test/PROPOSAL.cxx
)

target_link_libraries(PROPOSALtest PROPOSAL)

IF(DoTesting)
  execute_process(COMMAND mkdir ${PROPOSAL_BINARY_DIR}/bin/ OUTPUT_VARIABLE _output OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND  tar -xvf ${PROJECT_SOURCE_DIR}/Test/TestFiles.tar.Z -C ${PROPOSAL_BINARY_DIR}/bin/
                    OUTPUT_VARIABLE _output OUTPUT_STRIP_TRAILING_WHITESPACE)

  add_executable(UnitTest_StandardNormal Test/StandardNormal_Test.cxx)
  add_executable(UnitTest_Photonuclear Test/Photonuclear_TEST.cxx)
  add_executable(UnitTest_Integral Test/Integral_TEST.cxx)
  add_executable(UnitTest_Interpolant Test/Interpolant_TEST.cxx)
  add_executable(UnitTest_Bremsstrahlung Test/Bremsstrahlung_TEST.cxx)
  add_executable(UnitTest_Epairproduction Test/Epairproduction_TEST.cxx)
  add_executable(UnitTest_Ionization Test/Ionization_TEST.cxx)
  add_executable(UnitTest_RootFinder Test/RootFinder_TEST.cxx)
  add_executable(UnitTest_Medium Test/Medium_TEST.cxx)
  add_executable(UnitTest_Particle Test/Particle_TEST.cxx)
  add_executable(UnitTest_EnergyCutSettings Test/EnergyCutSettings_TEST.cxx)


  target_link_libraries(UnitTest_StandardNormal PROPOSAL gtest gtest_main)
  target_link_libraries(UnitTest_Integral PROPOSAL gtest gtest_main)
  target_link_libraries(UnitTest_Interpolant PROPOSAL gtest gtest_main)
  target_link_libraries(UnitTest_Bremsstrahlung PROPOSAL gtest gtest_main)
  target_link_libraries(UnitTest_Epairproduction PROPOSAL gtest gtest_main)
  target_link_libraries(UnitTest_Ionization PROPOSAL gtest gtest_main)
  target_link_libraries(UnitTest_Photonuclear PROPOSAL gtest gtest_main)
  target_link_libraries(UnitTest_RootFinder PROPOSAL gtest gtest_main)
  target_link_libraries(UnitTest_Medium PROPOSAL gtest gtest_main)
  target_link_libraries(UnitTest_Particle PROPOSAL gtest gtest_main)
  target_link_libraries(UnitTest_EnergyCutSettings PROPOSAL gtest gtest_main)

  add_test(UnitTest_RootFinder bin/UnitTest_RootFinder)
  add_test(UnitTest_StandardNormal bin/UnitTest_StandardNormal)
  add_test(UnitTest_Integral bin/UnitTest_Integral)
  add_test(UnitTest_Medium bin/UnitTest_Medium)
  add_test(UnitTest_Particle bin/UnitTest_Particle)
  add_test(UnitTest_EnergyCutSettings bin/UnitTest_EnergyCutSettings)
  add_test(UnitTest_Interpolant bin/UnitTest_Interpolant)
  add_test(UnitTest_Ionization bin/UnitTest_Ionization)
  add_test(UnitTest_Bremsstrahlung bin/UnitTest_Bremsstrahlung)
  add_test(UnitTest_Epairproduction bin/UnitTest_Epairproduction)
  add_test(UnitTest_Photonuclear bin/UnitTest_Photonuclear)

ENDIF()

ADD_SUBDIRECTORY( doc )
