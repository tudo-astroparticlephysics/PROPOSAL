cmake_minimum_required(VERSION 2.8)

project( PROPOSAL)


IF(I3_PROJECTS)
	add_definitions(-DICECUBE_PROJECT)
    include(cmake/I3_CMakeLists.cmake)
ELSE(I3_PROJECTS)
	add_definitions(-DPROPOSAL_STANDALONE)
    include(cmake/Standalone_CMakeLists.cmake)
ENDIF(I3_PROJECTS)

# Just for auto formating code
include(cmake/clang_cxx_dev_tools.cmake)
