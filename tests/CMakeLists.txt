add_subdirectory("${PROJECT_SOURCE_DIR}/vendor/google/googletest" "extern/googletest" EXCLUDE_FROM_ALL)
mark_as_advanced(
	BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
	gmock_build_tests gtest_build_samples gtest_build_tests
	gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
)


set_target_properties(gtest PROPERTIES FOLDER extern)
set_target_properties(gtest_main PROPERTIES FOLDER extern)
set_target_properties(gmock PROPERTIES FOLDER extern)
set_target_properties(gmock_main PROPERTIES FOLDER extern)

macro(package_add_test TESTNAME)
	add_executable(${TESTNAME} ${ARGN})
	target_link_libraries(${TESTNAME} gtest gmock gtest_main PROPOSAL)
	add_test(NAME ${TESTNAME} COMMAND ${TESTNAME} WORKING_DIRECTORY ${PROPOSAL_BINARY_DIR})
	set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()


ENABLE_TESTING()

EXECUTE_PROCESS(COMMAND mkdir -p ${PROPOSAL_BINARY_DIR}/bin/ OUTPUT_VARIABLE _output OUTPUT_STRIP_TRAILING_WHITESPACE)

#create tar directory with "tar -czvf TestFiles.tar.Z TestFiles/" and put it in Test directory
EXECUTE_PROCESS(COMMAND  tar -xvf ${PROJECT_SOURCE_DIR}/tests/TestFiles.tar.gz -C ${PROPOSAL_BINARY_DIR}/bin/
				OUTPUT_VARIABLE _output OUTPUT_STRIP_TRAILING_WHITESPACE)

package_add_test(UnitTest_Utility Utility_TEST.cxx)
package_add_test(UnitTest_Scattering Scattering_TEST.cxx)
package_add_test(UnitTest_Photonuclear Photonuclear_TEST.cxx)
package_add_test(UnitTest_Integral Integral_TEST.cxx)
package_add_test(UnitTest_Interpolant Interpolant_TEST.cxx)
package_add_test(UnitTest_Bremsstrahlung Bremsstrahlung_TEST.cxx)
package_add_test(UnitTest_Compton Compton_TEST.cxx)
package_add_test(UnitTest_Epairproduction Epairproduction_TEST.cxx)
package_add_test(UnitTest_Mupairproduction Mupairproduction_TEST.cxx)
package_add_test(UnitTest_WeakInteraction WeakInteraction_TEST.cxx)
package_add_test(UnitTest_Annihilation Annihilation_TEST.cxx)
package_add_test(UnitTest_Ionization Ionization_TEST.cxx)
package_add_test(UnitTest_PhotoPair PhotoPair_TEST.cxx)
package_add_test(UnitTest_Medium Medium_TEST.cxx)
package_add_test(UnitTest_Particle Particle_TEST.cxx)
package_add_test(UnitTest_ParticleDef ParticleDef_TEST.cxx)
package_add_test(UnitTest_DecayChannel DecayChannel_TEST.cxx)
package_add_test(UnitTest_DecayTable DecayTable_TEST.cxx)
package_add_test(UnitTest_EnergyCutSettings EnergyCutSettings_TEST.cxx)
package_add_test(UnitTest_ContinuousRandomization ContinuousRandomization_TEST.cxx)
package_add_test(UnitTest_Geometry Geometry_TEST.cxx)
package_add_test(UnitTest_Vector3D Vector3D_TEST.cxx)
package_add_test(UnitTest_Propagation Propagation_TEST.cxx)
package_add_test(UnitTest_Sector Sector_TEST.cxx)
package_add_test(UnitTest_MathMethods MathMethods_TEST.cxx)
package_add_test(UnitTest_Spline Spline_TEST.cxx)
package_add_test(UnitTest_Density Density_distribution_TEST.cxx)
