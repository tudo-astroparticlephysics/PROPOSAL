SET(CMAKE_CXX_STANDARD 17)

include(GoogleTest)

SET(RESOURCE_URL https://proposal.app.tu-dortmund.de/resources)
SET(TEST_FILES TestFiles-2021-01-13.tar.gz)

# download the test files
add_custom_command(OUTPUT TestFiles
    COMMAND curl --silent ${RESOURCE_URL}/${TEST_FILES} -o TestFiles.tar.gz
    COMMAND ${CMAKE_COMMAND} -E tar -xzf TestFiles.tar.gz
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Downloading and extracting test files."
    )

add_custom_target(testfiles DEPENDS TestFiles)

# small lib with utilities for PROPOSAL testing
add_library(TestUtilities INTERFACE)
add_library(PROPOSAL::TestUtilities ALIAS TestUtilities)
target_include_directories(TestUtilities INTERFACE PROPOSALTestUtilities)

find_package(GTest REQUIRED)
# creat test binary and link against proposal and gtest
macro(package_add_test TESTNAME)
    add_executable(${TESTNAME} ${ARGN})
    add_dependencies(${TESTNAME} testfiles)
    target_link_libraries(${TESTNAME} PROPOSAL::PROPOSAL PROPOSAL::TestUtilities GTest::GTest)
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
package_add_test(UnitTest_Secondaries Secondaries_TEST.cxx)
package_add_test(UnitTest_Sector Sector_TEST.cxx)
package_add_test(UnitTest_Utility Utility_TEST.cxx)
package_add_test(UnitTest_Time Time_TEST.cxx)
