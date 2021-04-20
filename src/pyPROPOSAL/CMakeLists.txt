cmake_minimum_required(VERSION 3.9)

project(pyPROPOSAL LANGUAGES CXX)

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
