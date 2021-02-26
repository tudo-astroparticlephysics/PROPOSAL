include(GNUInstallDirs)

add_library(PROPOSAL)
add_library(PROPOSAL::PROPOSAL ALIAS PROPOSAL)

target_include_directories(PROPOSAL PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )

add_subdirectory(PROPOSAL)
add_subdirectory(detail)

find_package(CubicInterpolation REQUIRED)
find_package(spdlog REQUIRED)

target_link_libraries(PROPOSAL PRIVATE
    CubicInterpolation::CubicInterpolation
    PUBLIC
    spdlog::spdlog
    )

install(TARGETS PROPOSAL EXPORT PROPOSALTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )

export(EXPORT PROPOSALTargets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/PROPOSALTargets.cmake")

install(EXPORT PROPOSALTargets
    FILE PROPOSALTargets.cmake
    NAMESPACE PROPOSAL::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PROPOSAL
    )

include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/PROPOSALConfig.cmake
    INSTALL_DESTINATION lib/cmake
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
    )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/PROPOSALConfig.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PROPOSAL
    )

file(GLOB_RECURSE Find_Modules "${PROJECT_BINARY_DIR}/Find*.cmake")

install(FILES ${Find_Modules}
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PROPOSAL)
