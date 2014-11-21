cmake_minimum_required(VERSION 2.8)

project( PROPOSAL)


IF(I3_PROJECTS)
    include(resources/I3_CMakeLists.cmake)
ELSE(I3_PROJECTS)
    include(resources/Standalone_CMakeLists.cmake)
ENDIF(I3_PROJECTS)
