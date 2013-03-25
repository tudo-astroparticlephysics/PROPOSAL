cmake_minimum_required(VERSION 2.6)

#FIND_LIBRARY(LOG4CPLUS_LIB NAMES log4cplus HINTS $LD_LIBRARY_PATH)
#MESSAGE(STATUS "log4cplus lib path: ${LOG4CPLUS_LIB}")
project( PROPOSAL)

IF(DEFINED ENV{GTEST_PATH})
  enable_testing()
  #1. Path
  #2. Path where the binarys will be copied to
  add_subdirectory($ENV{GTEST_PATH} ${CMAKE_CURRENT_BINARY_DIR}/gtestbin)
  set(DoTesting true)
  message("GTEST_PATH is set. Remember to add gtest/include to C_INCLUDE_PATH and CPP_INCLUDE_PATH")
ELSE()
  message("GTEST_PATH is not defined.  To build testsuite you have to set GTEST_PATH to gtest root directory")
ENDIF()



set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

configure_file(${PROJECT_SOURCE_DIR}/Test/Brems_dEdx.txt ${PROPOSAL_BINARY_DIR}/bin/Brems_dEdx.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/Test/Epair_dEdx.txt ${PROPOSAL_BINARY_DIR}/bin/Epair_dEdx.txt COPYONLY)


SET(DEBUG OFF)
option (DEBUG "DEBUG" OFF)
if(DEBUG)
        set( CMAKE_CXX_FLAGS "-g -pg -O0" )
endif()

include_directories("${PROJECT_SOURCE_DIR}/public")

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

  target_link_libraries(UnitTest_Integral PROPOSAL gtest gtest_main)
  target_link_libraries(UnitTest_Interpolant PROPOSAL gtest gtest_main)
  target_link_libraries(UnitTest_Bremsstrahlung PROPOSAL gtest gtest_main)
  target_link_libraries(UnitTest_Epairproduction PROPOSAL gtest gtest_main)

  add_test(UnitTest_Integral bin/UnitTest_Integral)
  add_test(UnitTest_Interpolant bin/UnitTest_Interpolant)
  add_test(UnitTest_Bremsstrahlung bin/UnitTest_Bremsstrahlung)
  add_test(UnitTest_Epairproduction bin/UnitTest_Epairproduction)

ENDIF()
