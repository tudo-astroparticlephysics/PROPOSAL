cmake_minimum_required(VERSION 2.6)

#FIND_LIBRARY(LOG4CPLUS_LIB NAMES log4cplus HINTS $LD_LIBRARY_PATH)
#MESSAGE(STATUS "log4cplus lib path: ${LOG4CPLUS_LIB}")

project( PROPOSAL)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories("${PROJECT_SOURCE_DIR}/public")

add_library(PROPOSAL
        private/PROPOSAL/Integral.cxx
        private/PROPOSAL/methods.cxx
        private/PROPOSAL/PhysicsModel.cxx
        private/PROPOSAL/MathModel.cxx
        private/PROPOSAL/Output.cxx
        private/PROPOSAL/Interpolant.cxx
        private/PROPOSAL/Bremsstrahlung.cxx
	)

set_target_properties(PROPOSAL PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror")



add_executable(PROPOSALtest
        private/test/PROPOSAL.cxx

)

target_link_libraries(PROPOSALtest PROPOSAL)
