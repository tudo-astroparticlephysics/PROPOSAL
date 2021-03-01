find_package(Doxygen REQUIRED)

set(DOXYGEN_GENERATE_HTML TRUE)
set(DOXYGEN_GENERATE_MAN FALSE)

file(GLOB_RECURSE PROPOSAL_HEADERS ${PROJECT_SOURCE_DIR}/src/PROPOSAL/*.h)

doxygen_add_docs(doxygen ALL
    ${PROPOSAL_HEADERS}
    COMMENT "Generate html documentation."
    )
