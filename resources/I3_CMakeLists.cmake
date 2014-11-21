# - Try to find Log4cplus
# Once done this will define
#  LOG4CPLUS_FOUND - System has Log4cplus
#  LOG4CPLUS_INCLUDE_DIRS - The Log4cplus include directories
#  LOG4CPLUS_LIBRARIES - The libraries needed to use Log4cplus
    add_definitions(-DPROPOSAL_STANDALONE=0)
i3_project(PROPOSAL)

i3_add_library(PROPOSAL
        private/PROPOSAL/Integral.cxx
        private/PROPOSAL/methods.cxx
        private/PROPOSAL/MathModel.cxx
        private/PROPOSAL/Bremsstrahlung.cxx
        private/PROPOSAL/CrossSections.cxx
        private/PROPOSAL/Decay.cxx
        private/PROPOSAL/Epairproduction.cxx
        private/PROPOSAL/Ionization.cxx
        private/PROPOSAL/Photonuclear.cxx
        private/PROPOSAL/Medium.cxx
        private/PROPOSAL/PROPOSALParticle.cxx
        private/PROPOSAL/EnergyCutSettings.cxx
        private/PROPOSAL/Interpolant.cxx
        private/PROPOSAL/StandardNormal.cxx
        private/PROPOSAL/RootFinder.cxx
        private/PROPOSAL/ProcessCollection.cxx
        private/PROPOSAL/Propagator.cxx
        private/PROPOSAL/ContinuousRandomization.cxx
        private/PROPOSAL/Geometry.cxx
        private/PROPOSAL/Scattering.cxx
        private/PROPOSAL/ScatteringFirstOrder.cxx
        private/PROPOSAL/ScatteringMoliere.cxx
        private/PROPOSAL/Output.cxx
        private/PROPOSAL-icetray/I3PropagatorServicePROPOSAL.cxx
        private/PROPOSAL-icetray/SimplePropagator.cxx

        USE_TOOLS boost
        USE_PROJECTS icetray dataclasses sim-services simclasses phys-services
        )

set_target_properties(PROPOSAL PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS}  -Wall")

i3_add_pybindings(PROPOSAL
        private/PROPOSAL-icetray/pybindings.cxx
        USE_PROJECTS PROPOSAL
)

# Pre-generate parameterization tables for use in read-only environments
add_custom_command(OUTPUT resources/tables/.tables.auto_generated COMMAND ${CMAKE_BINARY_DIR}/env-shell.sh python resources/tables/generate_tables.py DEPENDS icetray-pybindings sim_services-pybindings PROPOSAL-pybindings)
add_custom_target(PROPOSAL-tables ALL DEPENDS resources/tables/.tables.auto_generated)

set(LIB_${PROJECT_NAME}_TESTS
        private/PROPOSAL-icetray/test/main.cxx
)

if (SPRNG_FOUND)
        # this test requires SPRNG
        LIST(APPEND LIB_${PROJECT_NAME}_TESTS
                private/PROPOSAL-icetray/test/Repeatability.cxx
        )
endif (SPRNG_FOUND)


i3_test_executable(test
        ${LIB_${PROJECT_NAME}_TESTS}
        USE_TOOLS boost gsl
        USE_PROJECTS PROPOSAL icetray dataclasses phys-services
)

