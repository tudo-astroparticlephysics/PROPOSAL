i3_project( PROPOSAL)

i3_add_library(PROPOSAL
	private/PROPOSAL/Integral.cxx
	private/PROPOSAL/DFunctionOfx.cxx
	private/PROPOSAL/Interpolate.cxx
	private/PROPOSAL/MathModel.cxx
	private/PROPOSAL/FunctionOfx.cxx
	private/PROPOSAL/PhysicsModel.cxx
	private/PROPOSAL/Medium.cxx
	private/PROPOSAL/methods.cxx
	private/PROPOSAL/Output.cxx
	private/PROPOSAL/Ionizationloss.cxx
	private/PROPOSAL/IonizStochastic.cxx
	private/PROPOSAL/IonizContinuous.cxx
	private/PROPOSAL/Epairproduction.cxx
	private/PROPOSAL/EpairContinuous.cxx
	private/PROPOSAL/EpairStochastic.cxx
	private/PROPOSAL/Photonuclear.cxx
	private/PROPOSAL/PhotoStochastic.cxx
	private/PROPOSAL/PhotoContinuous.cxx
	private/PROPOSAL/Bremsstrahlung.cxx
	private/PROPOSAL/BremsStochastic.cxx
	private/PROPOSAL/BremsContinuous.cxx
	private/PROPOSAL/Decay.cxx
	private/PROPOSAL/FindRoot.cxx
	private/PROPOSAL/Energy2Loss.cxx
	private/PROPOSAL/Energy2LossE.cxx
	private/PROPOSAL/Energy2LossX.cxx
	private/PROPOSAL/StandardNormal.cxx
	private/PROPOSAL/Scattering.cxx
	private/PROPOSAL/CrossSections.cxx
	private/PROPOSAL/PROPOSALParticle.cxx
	private/PROPOSAL/Propagate.cxx
	private/PROPOSAL/MMCexception.cxx
	private/PROPOSAL/Mpairproduction.cxx
	private/PROPOSAL/MpairContinuous.cxx
	private/PROPOSAL/MpairStochastic.cxx
	private/PROPOSAL/Amanda.cxx
	private/PROPOSAL/FunctionInt2.cxx

	private/PROPOSAL-icetray/I3PropagatorServicePROPOSAL.cxx
	private/PROPOSAL-icetray/I3PropagatorPROPOSAL.cxx

	USE_TOOLS boost
	USE_PROJECTS icetray dataclasses sim-services simclasses   
	)

i3_add_pybindings(PROPOSAL
        private/PROPOSAL-icetray/pybindings.cxx
        USE_PROJECTS PROPOSAL
)

i3_test_executable(PROPOSALtest 
	private/test/PROPOSAL.cxx
	USE_PROJECTS PROPOSAL
)

