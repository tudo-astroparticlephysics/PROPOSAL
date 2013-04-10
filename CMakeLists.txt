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

configure_file(${PROJECT_SOURCE_DIR}/Test/Brems_dEdx.txt ${PROPOSAL_BINARY_DIR}/bin/Brems_dEdx.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/Test/Brems_dNdx.txt ${PROPOSAL_BINARY_DIR}/bin/Brems_dNdx.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/Test/Brems_dNdxrnd.txt ${PROPOSAL_BINARY_DIR}/bin/Brems_dNdxrnd.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/Test/Brems_e.txt ${PROPOSAL_BINARY_DIR}/bin/Brems_e.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/Test/Brems_dNdxrnd_interpol.txt ${PROPOSAL_BINARY_DIR}/bin/Brems_dNdxrnd_interpol.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/Test/Brems_dEdx_interpol.txt ${PROPOSAL_BINARY_DIR}/bin/Brems_dEdx_interpol.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/Test/Brems_dNdx_interpol.txt ${PROPOSAL_BINARY_DIR}/bin/Brems_dNdx_interpol.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/Test/Brems_e_interpol.txt ${PROPOSAL_BINARY_DIR}/bin/Brems_e_interpol.txt COPYONLY)

configure_file(${PROJECT_SOURCE_DIR}/Test/Epair_dEdx.txt ${PROPOSAL_BINARY_DIR}/bin/Epair_dEdx.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/Test/Epair_dNdx.txt ${PROPOSAL_BINARY_DIR}/bin/Epair_dNdx.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/Test/Epair_dNdxrnd.txt ${PROPOSAL_BINARY_DIR}/bin/Epair_dNdxrnd.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/Test/Epair_e.txt ${PROPOSAL_BINARY_DIR}/bin/Epair_e.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/Test/Epair_dEdx_interpol.txt ${PROPOSAL_BINARY_DIR}/bin/Epair_dEdx_interpol.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/Test/Epair_dNdx_interpol.txt ${PROPOSAL_BINARY_DIR}/bin/Epair_dNdx_interpol.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/Test/Epair_dNdxrnd_interpol.txt ${PROPOSAL_BINARY_DIR}/bin/Epair_dNdxrnd_interpol.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/Test/Epair_e_interpol.txt ${PROPOSAL_BINARY_DIR}/bin/Epair_e_interpol.txt COPYONLY)

configure_file(${PROJECT_SOURCE_DIR}/Test/Ioniz_dEdx.txt ${PROPOSAL_BINARY_DIR}/bin/Ioniz_dEdx.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/Test/Ioniz_dNdx.txt ${PROPOSAL_BINARY_DIR}/bin/Ioniz_dNdx.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/Test/Ioniz_dNdxrnd.txt ${PROPOSAL_BINARY_DIR}/bin/Ioniz_dNdxrnd.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/Test/Ioniz_e.txt ${PROPOSAL_BINARY_DIR}/bin/Ioniz_e.txt COPYONLY)

configure_file(${PROJECT_SOURCE_DIR}/Test/Photo_dEdx.txt ${PROPOSAL_BINARY_DIR}/bin/Photo_dEdx.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/Test/Photo_dNdx.txt ${PROPOSAL_BINARY_DIR}/bin/Photo_dNdx.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/Test/Photo_dNdxrnd.txt ${PROPOSAL_BINARY_DIR}/bin/Photo_dNdxrnd.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/Test/Photo_e.txt ${PROPOSAL_BINARY_DIR}/bin/Photo_e.txt COPYONLY)

configure_file(${PROJECT_SOURCE_DIR}/Test/rnd.txt ${PROPOSAL_BINARY_DIR}/bin/rnd.txt COPYONLY)

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
	)

set_target_properties(PROPOSAL PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror")

add_executable(PROPOSALtest
        private/test/PROPOSAL.cxx
)

target_link_libraries(PROPOSALtest PROPOSAL)

IF(DoTesting)
  add_executable(UnitTest_Integral Test/Integral_TEST.cxx)
  add_executable(UnitTest_Interpolant Test/Interpolant_TEST.cxx)
  add_executable(UnitTest_Bremsstrahlung Test/Bremsstrahlung_TEST.cxx)
  add_executable(UnitTest_Epairproduction Test/Epairproduction_TEST.cxx)
  add_executable(UnitTest_Ionization Test/Ionization_TEST.cxx)
  add_executable(UnitTest_Photonuclear Test/Photonuclear_TEST.cxx)

  target_link_libraries(UnitTest_Integral PROPOSAL gtest gtest_main)
  target_link_libraries(UnitTest_Interpolant PROPOSAL gtest gtest_main)
  target_link_libraries(UnitTest_Bremsstrahlung PROPOSAL gtest gtest_main)
  target_link_libraries(UnitTest_Epairproduction PROPOSAL gtest gtest_main)
  target_link_libraries(UnitTest_Ionization PROPOSAL gtest gtest_main)
  target_link_libraries(UnitTest_Photonuclear PROPOSAL gtest gtest_main)

  add_test(UnitTest_Integral bin/UnitTest_Integral)
  add_test(UnitTest_Interpolant bin/UnitTest_Interpolant)
  add_test(UnitTest_Bremsstrahlung bin/UnitTest_Bremsstrahlung)
  add_test(UnitTest_Epairproduction bin/UnitTest_Epairproduction)
  add_test(UnitTest_Ionization bin/UnitTest_Ionization)
  add_test(UnitTest_Photonuclear bin/UnitTest_Photonuclear)

ENDIF()
