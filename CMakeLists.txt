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
ELSE()
  message("GTEST_PATH is not defined.  To build testsuite you have to set GTEST_PATH to gtest root directory")
ENDIF()



set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories("${PROJECT_SOURCE_DIR}/public")

add_library(PROPOSAL
        private/PROPOSAL/Integral.cxx
        private/PROPOSAL/methods.cxx
        private/PROPOSAL/MathModel.cxx
        private/PROPOSAL/Output.cxx
        #private/PROPOSAL/Interpolant.cxx
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
  add_executable(UnitTest
	  Test/Integral_TEST.cxx
  )
  target_link_libraries(UnitTest PROPOSAL gtest gtest_main)
  add_test(UnitTests bin/UnitTest)
ENDIF()