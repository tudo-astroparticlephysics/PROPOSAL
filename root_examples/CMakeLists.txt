# Load some basic macros which are needed later on
include(FindROOT.cmake)

IF(ROOT_FOUND)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/root_examples)

    set(INCLUDE_DIRECTORIES
    ${ROOT_INCLUDE_DIR}
    )

    include_directories( ${INCLUDE_DIRECTORIES})

    set(LINK_DIRECTORIES
    ${ROOT_LIBRARY_DIR}
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    )

    link_directories( ${LINK_DIRECTORIES})

    add_executable(Plot_range plot_range.cxx)
    add_executable(Plot_Ionization_dEdx plot_Ionization_dEdx.cxx)
    add_executable(Plot_Ionization_dNdx plot_Ionization_dNdx.cxx)
    add_executable(Plot_Epairproduction_dEdx plot_Epairproduction_dEdx.cxx)

    target_link_libraries(Plot_range ${ROOT_LIBRARIES})
    target_link_libraries(Plot_range PROPOSAL)

    target_link_libraries(Plot_Ionization_dEdx ${ROOT_LIBRARIES})
    target_link_libraries(Plot_Ionization_dEdx PROPOSAL)

    target_link_libraries(Plot_Ionization_dNdx ${ROOT_LIBRARIES})
    target_link_libraries(Plot_Ionization_dNdx PROPOSAL)

    target_link_libraries(Plot_Epairproduction_dEdx ${ROOT_LIBRARIES})
    target_link_libraries(Plot_Epairproduction_dEdx PROPOSAL)

ELSE(ROOT_FOUND)
    MESSAGE(STATUS "ROOT not found...")
    MESSAGE(STATUS "ROOT examples will not be builded.")
    MESSAGE(STATUS "Make sure you have ROOT installed and ROOTSYS is set.")

ENDIF(ROOT_FOUND)

