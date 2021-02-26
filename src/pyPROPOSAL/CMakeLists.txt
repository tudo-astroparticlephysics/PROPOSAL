message(STATUS "Building the python wrapper library.")
find_package(PythonLibs REQUIRED)
find_package(pybind11 REQUIRED)
file(GLOB_RECURSE PYTHON_SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cxx)

pybind11_add_module(pyproposal SHARED ${PYTHON_SRC_FILES})
set_target_properties(pyproposal PROPERTIES OUTPUT_NAME proposal)
target_include_directories(pyproposal PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:${PYTHON_INCLUDE_DIRS}>
    )

target_link_libraries(pyproposal PRIVATE PROPOSAL)
target_compile_options(pyproposal PRIVATE -fvisibility=hidden)
install(TARGETS pyproposal EXPORT PROPOSALTargets
    DESTINATION ${CMAKE_INSTALL_LIBDIR})
