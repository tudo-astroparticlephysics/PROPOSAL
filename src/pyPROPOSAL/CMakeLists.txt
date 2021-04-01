cmake_minimum_required(VERSION 3.9)

project(pyPROPOSAL LANGUAGES CXX)

if (EXISTS ${CMAKE_BINARY_DIR}/Findpybind11.cmake)
    message(WARNING " There seems to be currently a bug in conans pybind11 module, you have to temporary install pybind11 by your own and delete the Findpybind11.cmake module in your binary dir. Hopefully fix this soon.")
endif()

find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)


file(GLOB_RECURSE PYTHON_SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cxx)
pybind11_add_module(pyPROPOSAL SHARED
    ${PYTHON_SRC_FILES})

set_target_properties(pyPROPOSAL PROPERTIES
    OUTPUT_NAME proposal)

target_include_directories(pyPROPOSAL PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${PYTHON_INCLUDE_DIRS}>
    )

target_link_libraries(pyPROPOSAL PRIVATE PROPOSAL::PROPOSAL)
install(TARGETS pyPROPOSAL EXPORT pyPROPOSALTargets
    DESTINATION ${CMAKE_INSTALL_LIBDIR})
