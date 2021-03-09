cmake_minimum_required(VERSION 3.8)

project(PROPOSAL VERSION 6.1.6 LANGUAGES CXX)


IF(I3_PROJECTS)
	add_definitions(-DICECUBE_PROJECT)
    include(cmake/I3_CMakeLists.cmake)
ELSE(I3_PROJECTS)
	add_definitions(-DPROPOSAL_STANDALONE)
    include(cmake/Standalone_CMakeLists.cmake)
ENDIF(I3_PROJECTS)

# Just for auto formating code
include(cmake/clang_cxx_dev_tools.cmake)
