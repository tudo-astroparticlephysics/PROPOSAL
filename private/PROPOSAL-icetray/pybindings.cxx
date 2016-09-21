
#include "PROPOSAL-icetray/I3PropagatorServicePROPOSAL.h"
#include "icetray/load_project.h"

I3_PYTHON_MODULE(PROPOSAL)
{
	load_project("PROPOSAL", false);
	
	using namespace boost::python;
	import("icecube.dataclasses");
	import("icecube.sim_services");
	
	enum_<I3PropagatorServicePROPOSAL::BremsstrahlungParametrization>("BremsstrahlungParametrization")
		.value("KelnerKokoulinPetrukhin", I3PropagatorServicePROPOSAL::KelnerKokoulinPetrukhin)
		.value("PetrukhinShestakov", I3PropagatorServicePROPOSAL::PetrukhinShestakov)
		.value("AndreevBerzrukovBugaev", I3PropagatorServicePROPOSAL::AndreevBerzrukovBugaev)
		.value("CompleteScreeningCase", I3PropagatorServicePROPOSAL::CompleteScreeningCase)
	;
	
	enum_<I3PropagatorServicePROPOSAL::PhotonuclearParametrizationFamily>("PhotonuclearParametrizationFamily")
		.value("BezrukovBugaevSoft", I3PropagatorServicePROPOSAL::BezrukovBugaevSoft)
		.value("BezrukovBugaevHard", I3PropagatorServicePROPOSAL::BezrukovBugaevHard)
		.value("AbramowiczLevinLevyMaor", I3PropagatorServicePROPOSAL::AbramowiczLevinLevyMaor)
		.value("ButkevichMikheyev", I3PropagatorServicePROPOSAL::ButkevichMikheyev)
	;
	
	enum_<I3PropagatorServicePROPOSAL::PhotonuclearParametrization>("PhotonuclearParametrization")
		.value("AbramowiczLevinLevyMaor91", I3PropagatorServicePROPOSAL::AbramowiczLevinLevyMaor91)
		.value("AbramowiczLevinLevyMaor97", I3PropagatorServicePROPOSAL::AbramowiczLevinLevyMaor97)
		.value("BezrukovBugaev", I3PropagatorServicePROPOSAL::BezrukovBugaev)
		.value("ZEUS", I3PropagatorServicePROPOSAL::ZEUS)
	;
	
	enum_<I3PropagatorServicePROPOSAL::ShadowingParametrization>("ShadowingParametrization")
		.value("Dutta", I3PropagatorServicePROPOSAL::Dutta)
		.value("Butkevich", I3PropagatorServicePROPOSAL::Butkevich)
	;
		
	class_<I3PropagatorServicePROPOSAL, boost::shared_ptr<I3PropagatorServicePROPOSAL>,
	    bases<I3PropagatorService>, boost::noncopyable>(
	    "I3PropagatorServicePROPOSAL",
	    init<std::string,std::string,double,double,I3Particle::ParticleType,double,
	    I3PropagatorServicePROPOSAL::BremsstrahlungParametrization,
	    I3PropagatorServicePROPOSAL::PhotonuclearParametrizationFamily,
	    I3PropagatorServicePROPOSAL::PhotonuclearParametrization,
	    I3PropagatorServicePROPOSAL::ShadowingParametrization>(
	    (arg("mediadef")=I3PropagatorServicePROPOSAL::GetDefaultMediaDef(),
	     arg("tabledir")=I3PropagatorServicePROPOSAL::GetDefaultTableDir(),
	     arg("cylinderRadius")=800*I3Units::m, arg("cylinderHeight")=1600*I3Units::m,
	     arg("type")=I3Particle::MuMinus, arg("particleMass")=NAN,
	     arg("bremsstrahlungParametrization")=I3PropagatorServicePROPOSAL::KelnerKokoulinPetrukhin,
	     arg("photonuclearParametrizationFamily")=I3PropagatorServicePROPOSAL::AbramowiczLevinLevyMaor,
	     arg("photonuclearParametrization")=I3PropagatorServicePROPOSAL::AbramowiczLevinLevyMaor97,
	     arg("nuclearShadowingParametrization")=I3PropagatorServicePROPOSAL::Butkevich), 
	     ":param mediadef: Path the the media definition file\n"
	     ":param tabledir: Path to a directory in which to store interpolation constants for cross-section integrals\n"
	     ":param cylinderRadius: Radius of the target volume in meters\n"
	     ":param cylinderHeight: Full height of the target volume in meters\n"
	     ":param type: Type of particle to propagate\n"
	     ":param particleMass: Mass of the propagated particle in GeV. This is only used if type is something exotic.\n"
	     ":param bremsstrahlungParametrization: Parametrization of the bremsstrahlung cross-section to use\n"
	     ":param photonuclearParametrizationFamily: Family of photonuclear cross-section parametrization to use\n"
	     ":param photonuclearParametrization: Specific edition of the photonuclear cross-section parametrization to use\n"
	     ":param nuclearShadowingParametrization: Nuclear shadowing parametrization to use\n"))
	    .def("set_tear_down_per_call", &I3PropagatorServicePROPOSAL::SetTearDownPerCall)
	;

}
