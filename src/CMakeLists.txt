include(GNUInstallDirs)

if(SKBUILD)
    message(STATUS "SKBUILD")
    set(CMAKE_INSTALL_LIBDIR ".")
    set(CMAKE_INSTALL_INCLUDEDIR "proposal/include")
endif()
add_subdirectory(PROPOSAL)

if(BUILD_PYTHON)
    add_subdirectory(pyPROPOSAL)
endif()
