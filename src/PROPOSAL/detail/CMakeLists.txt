file(GLOB_RECURSE PROPOSAL_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cxx)

target_sources(PROPOSAL PRIVATE
    ${PROPOSAL_sources})
