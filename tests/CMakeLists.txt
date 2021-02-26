include(GoogleTest)
enable_testing()

include(FetchContent)
file(DOWNLOAD
    https://proposal.app.tu-dortmund.de/resources/TestFiles_2020-10-21.tar.gz
    ${CMAKE_CURRENT_BINARY_DIR}/TestFiles.tar.gz)
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xvzf TestFiles.tar.gz
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

find_package(GTest REQUIRED)

macro(package_add_test TESTNAME)
    add_executable(${TESTNAME} ${ARGN})
    target_link_libraries(${TESTNAME} PROPOSAL::PROPOSAL GTest::GTest)
    gtest_discover_tests(${TESTNAME})
endmacro()

# Unit Tests
package_add_test(UnitTest_Decay Decay_TEST.cxx)
package_add_test(UnitTest_DecayChannel DecayChannel_TEST.cxx)
package_add_test(UnitTest_DecayTable DecayTable_TEST.cxx)
package_add_test(UnitTest_Density Density_distribution_TEST.cxx)
package_add_test(UnitTest_EnergyCutSettings EnergyCutSettings_TEST.cxx)
package_add_test(UnitTest_Geometry Geometry_TEST.cxx)
package_add_test(UnitTest_Integral Integral_TEST.cxx)
package_add_test(UnitTest_Interpolant Interpolant_TEST.cxx)
package_add_test(UnitTest_MathMethods MathMethods_TEST.cxx)
package_add_test(UnitTest_Medium Medium_TEST.cxx)
package_add_test(UnitTest_Particle Particle_TEST.cxx)
package_add_test(UnitTest_ParticleDef ParticleDef_TEST.cxx)
package_add_test(UnitTest_Spline Spline_TEST.cxx)
package_add_test(UnitTest_Vector3D Vector3D_TEST.cxx)

# cross section tests
package_add_test(UnitTest_Annihilation Annihilation_TEST.cxx)
package_add_test(UnitTest_Bremsstrahlung Bremsstrahlung_TEST.cxx)
package_add_test(UnitTest_Compton Compton_TEST.cxx)
package_add_test(UnitTest_Epairproduction Epairproduction_TEST.cxx)
package_add_test(UnitTest_Ionization Ionization_TEST.cxx)
package_add_test(UnitTest_Mupairproduction Mupairproduction_TEST.cxx)
package_add_test(UnitTest_Photonuclear Photonuclear_TEST.cxx)
package_add_test(UnitTest_PhotoPair PhotoPair_TEST.cxx)
package_add_test(UnitTest_WeakInteraction WeakInteraction_TEST.cxx)

# propagation and utility tests
package_add_test(UnitTest_ContinuousRandomization ContinuousRandomization_TEST.cxx)
package_add_test(UnitTest_Displacement Displacement_TEST.cxx)
package_add_test(UnitTest_Interaction Interaction_TEST.cxx)
package_add_test(UnitTest_PropagationUtilityIntegral UtilityIntegral_TEST.cxx)
package_add_test(UnitTest_PropagationUtilityInterpolant UtilityInterpolant_TEST.cxx)
package_add_test(UnitTest_Scattering Scattering_TEST.cxx)
# package_add_test(UnitTest_Secondaries Secondaries_TEST.cxx)
package_add_test(UnitTest_Sector Sector_TEST.cxx)
package_add_test(UnitTest_Utility Utility_TEST.cxx)
package_add_test(UnitTest_Time Time_TEST.cxx)
