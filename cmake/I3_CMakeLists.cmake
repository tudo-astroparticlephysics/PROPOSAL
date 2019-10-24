
i3_project(PROPOSAL
    DOCS_DIR resources/doc
)

# file(GLOB_RECURSE PROPOSAL_SRC_FILES ${PROJECT_SOURCE_DIR}/private/PROPOSAL/*)
set (PROPOSAL_SRC_FILES
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/Constants.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/EnergyCutSettings.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/Output.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/Propagator.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/PropagatorService.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/ComptonIntegral.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/ComptonInterpolant.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/BremsIntegral.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/BremsInterpolant.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/CrossSection.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/CrossSectionIntegral.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/CrossSectionInterpolant.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/EpairIntegral.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/EpairInterpolant.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/MupairIntegral.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/MupairInterpolant.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/WeakIntegral.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/WeakInterpolant.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/PhotoPairIntegral.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/PhotoPairInterpolant.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/IonizIntegral.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/IonizInterpolant.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/PhotoIntegral.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/PhotoInterpolant.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/factories/BremsstrahlungFactory.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/factories/ComptonFactory.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/factories/EpairProductionFactory.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/factories/PhotoPairFactory.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/factories/MupairProductionFactory.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/factories/IonizationFactory.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/factories/PhotonuclearFactory.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/factories/WeakInteractionFactory.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/parametrization/Bremsstrahlung.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/parametrization/Compton.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/parametrization/EpairProduction.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/parametrization/MupairProduction.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/parametrization/PhotoPairProduction.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/parametrization/WeakInteraction.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/parametrization/ParamTables.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/parametrization/Ionization.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/parametrization/Parametrization.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/parametrization/PhotoQ2Integration.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/parametrization/PhotoRealPhotonAssumption.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/crossection/parametrization/Photonuclear.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/decay/DecayChannel.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/decay/DecayTable.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/decay/LeptonicDecayChannel.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/decay/ManyBodyPhaseSpace.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/decay/StableChannel.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/decay/TwoBodyPhaseSpace.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/geometry/Box.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/geometry/Cylinder.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/geometry/Geometry.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/geometry/GeometryFactory.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/geometry/Sphere.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/math/Integral.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/math/Interpolant.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/math/MathMethods.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/math/InterpolantBuilder.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/math/RandomGenerator.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/math/Vector3D.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/medium/Components.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/medium/Medium.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/medium/MediumFactory.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/methods.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/particle/Particle.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/particle/ParticleDef.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/propagation_utility/ContinuousRandomizer.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/propagation_utility/PropagationUtility.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/propagation_utility/PropagationUtilityIntegral.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/propagation_utility/PropagationUtilityInterpolant.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/scattering/Coefficients.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/scattering/Scattering.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/scattering/ScatteringFactory.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/scattering/ScatteringHighland.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/scattering/ScatteringHighlandIntegral.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/scattering/ScatteringMoliere.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/scattering/ScatteringNoScattering.cxx
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL/sector/Sector.cxx
)


i3_add_library(PROPOSAL
    ${PROPOSAL_SRC_FILES}
    private/PROPOSAL-icetray/I3PropagatorServicePROPOSAL.cxx
    private/PROPOSAL-icetray/SimplePropagator.cxx
    private/PROPOSAL-icetray/Converter.cxx

    USE_TOOLS boost
    USE_PROJECTS icetray serialization dataclasses sim-services simclasses phys-services
)

set_target_properties(PROPOSAL PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS}  -Wall -Werror -std=c++11")

i3_add_pybindings(PROPOSAL
    private/PROPOSAL-icetray/pybindings.cxx
    USE_TOOLS boost python
    USE_PROJECTS PROPOSAL
)

# Pre-generate parameterization tables for use in read-only environments
add_executable(PROPOSAL_table_creation
    ${PROJECT_SOURCE_DIR}/private/PROPOSAL-icetray/table_creation.cxx
)
set_target_properties(PROPOSAL_table_creation PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS}  -Wall -Werror -std=c++11")
target_link_libraries(PROPOSAL_table_creation PROPOSAL)

add_custom_command(
    OUTPUT ${PROJECT_SOURCE_DIR}/resources/tables/.tables.auto_generated
    COMMAND ${CMAKE_BINARY_DIR}/env-shell.sh ${CMAKE_BINARY_DIR}/bin/PROPOSAL_table_creation
    DEPENDS icetray PROPOSAL PROPOSAL_table_creation
)
add_custom_target(
    PROPOSAL_tables 
    DEPENDS ${CMAKE_BINARY_DIR}/bin/PROPOSAL_table_creation
    ${PROJECT_SOURCE_DIR}/resources/tables/.tables.auto_generated
)
add_dependencies(PROPOSAL_tables PROPOSAL_table_creation)

add_custom_command(TARGET PROPOSAL
    POST_BUILD
    COMMENT "To create PROPOSAL tables run: make PROPOSAL_tables"
    COMMAND echo "***************************************************************************"
    COMMAND echo "***                                                                     ***"
    COMMAND echo "***                  To build tables for PROPOSAL run:                  ***"
    COMMAND echo "***                        make PROPOSAL_tables                         ***"
    COMMAND echo "***         enter the ICETRAY environment and run the command:          ***"
    COMMAND echo "***                $I3_BUILD/bin/PROPOSAL_table_creation                ***"
    COMMAND echo "***                                                                     ***"
    COMMAND echo "***************************************************************************"
    USES_TERMINAL VERBATIM
    )

set(LIB_${PROJECT_NAME}_TESTS
    private/PROPOSAL-icetray/test/main.cxx
)

# FIXME: See https://code.icecube.wisc.edu/projects/icecube/ticket/2194
if (SPRNG_FOUND)
   # this test requires SPRNG
   LIST(APPEND LIB_${PROJECT_NAME}_TESTS
       private/PROPOSAL-icetray/test/Repeatability.cxx
   )
endif (SPRNG_FOUND)


i3_test_executable(test
    ${LIB_${PROJECT_NAME}_TESTS}
    USE_TOOLS boost gsl python
    USE_PROJECTS PROPOSAL icetray dataclasses phys-services
)
