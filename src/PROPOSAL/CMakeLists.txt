include(GNUInstallDirs)

add_library(PROPOSAL)
add_library(PROPOSAL::PROPOSAL ALIAS PROPOSAL)

if(MSVC)
    target_compile_options(PROPOSAL PRIVATE "/bigobj" "/EHsc")
elseif(MINGW)
    target_compile_options(PROPOSAL PRIVATE "-Wa,-mbig-obj")
endif()

target_include_directories(PROPOSAL PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )

find_package(CubicInterpolation REQUIRED)
find_package(spdlog REQUIRED)
find_package(nlohmann_json REQUIRED)

add_subdirectory(PROPOSAL)
add_subdirectory(detail)

target_link_libraries(PROPOSAL
    CubicInterpolation::CubicInterpolation
    spdlog::spdlog
    nlohmann_json::nlohmann_json
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

file(GLOB CMAKE_FILES "${CMAKE_BINARY_DIR}/*.cmake")
install(FILES ${CMAKE_FILES} DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PROPOSAL)
