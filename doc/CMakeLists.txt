# add a target to generate PROPOSAL documentation with Doxygen
SET( DOXYGEN_SOURCE_DIR
${CMAKE_CURRENT_SOURCE_DIR}/../src
)

find_package(Doxygen)
if(DOXYGEN_FOUND)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/images/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/images)

set(DOC_IMAGES ${CMAKE_CURRENT_BINARY_DIR}/images)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile )
add_custom_target(doc
${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
COMMENT "Generating PROPOSAL documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)


