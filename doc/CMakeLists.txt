# add a target to generate PROPOSAL documentation with Doxygen
SET( DOXYGEN_SOURCE_DIR
${CMAKE_CURRENT_SOURCE_DIR}/../src
)

find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc
${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
COMMENT "Generating PROPOSAL documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)


