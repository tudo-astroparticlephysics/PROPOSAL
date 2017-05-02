
#include "PROPOSAL-icetray/I3PropagatorServicePROPOSAL.h"
#include "PROPOSAL-icetray/SimplePropagator.h"
#include "icetray/load_project.h"

I3_PYTHON_MODULE(PROPOSAL)
{
	load_project("PROPOSAL", false);

	using namespace boost::python;
	import("icecube.dataclasses");
	import("icecube.sim_services");

    enum_<ParametrizationType::Enum>("BremsstrahlungParametrization")
        .value("KelnerKokoulinPetrukhin" , ParametrizationType::BremsKelnerKokoulinPetrukhin)
        .value("PetrukhinShestakov"      ,  ParametrizationType::BremsPetrukhinShestakov)
        .value("AndreevBezrrukovBugaev"  ,  ParametrizationType::BremsAndreevBezrukovBugaev)
        .value("CompleteScreeningCase"   ,  ParametrizationType::BremsCompleteScreeningCase)
    ;

    enum_<ParametrizationType::Enum>("PhotonuclearParametrization")
        .value("KokoulinShadowBezrukovSoft"               , ParametrizationType::PhotoKokoulinShadowBezrukovSoft)
        .value("KokoulinShadowBezrukovHard"               , ParametrizationType::PhotoKokoulinShadowBezrukovHard)
        .value("RhodeShadowBezrukovSoft"                  , ParametrizationType::PhotoRhodeShadowBezrukovSoft)
        .value("RhodeShadowBezrukovHard"                  , ParametrizationType::PhotoRhodeShadowBezrukovHard)
        .value("BezrukovBugaevShadowBezrukovSoft"         , ParametrizationType::PhotoBezrukovBugaevShadowBezrukovSoft)
        .value("BezrukovBugaevShadowBezrukovHard"         , ParametrizationType::PhotoBezrukovBugaevShadowBezrukovHard)
        .value("ZeusShadowBezrukovSoft"                   , ParametrizationType::PhotoZeusShadowBezrukovSoft)
        .value("ZeusShadowBezrukovHard"                   , ParametrizationType::PhotoZeusShadowBezrukovHard)
        .value("AbramowiczLevinLevyMaor91ShadowDutta"     , ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowDutta)
        .value("AbramowiczLevinLevyMaor91ShadowButkevich" , ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowButkevich)
        .value("AbramowiczLevinLevyMaor97ShadowDutta"     , ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowDutta)
        .value("AbramowiczLevinLevyMaor97ShadowButkevich" , ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich)
        .value("ButkevichMikhailovShadowDutta"            , ParametrizationType::PhotoButkevichMikhailovShadowDutta)
        .value("ButkevichMikhailovShadowButkevich"        , ParametrizationType::PhotoButkevichMikhailovShadowButkevich)
    ;


	class_<I3PropagatorServicePROPOSAL, boost::shared_ptr<I3PropagatorServicePROPOSAL>,
	    bases<I3PropagatorService>, boost::noncopyable>(
	    "I3PropagatorServicePROPOSAL",
	    init<std::string,std::string,double,double,I3Particle::ParticleType,double,
	    ParametrizationType::Enum,
	    ParametrizationType::Enum>(
	    (arg("mediadef")=I3PropagatorServicePROPOSAL::GetDefaultMediaDef(),
	     arg("tabledir")=I3PropagatorServicePROPOSAL::GetDefaultTableDir(),
	     arg("cylinderRadius")=800*I3Units::m,
	     arg("cylinderHeight")=1600*I3Units::m,
	     arg("type")=I3Particle::MuMinus,
	     arg("particleMass")=NAN,
	     arg("bremsstrahlungParametrization")=ParametrizationType::BremsKelnerKokoulinPetrukhin,
	     arg("photonuclearParametrization")=ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich),
	     ":param mediadef: Path the the media definition file\n"
	     ":param tabledir: Path to a directory in which to store interpolation constants for cross-section integrals\n"
	     ":param cylinderRadius: Radius of the target volume in meters\n"
	     ":param cylinderHeight: Full height of the target volume in meters\n"
	     ":param type: Type of particle to propagate\n"
	     ":param particleMass: Mass of the propagated particle in GeV. This is only used if type is something exotic.\n"
	     ":param bremsstrahlungParametrization: Parametrization of the bremsstrahlung cross-section to use\n"
	     ":param photonuclearParametrization: Parametrization of photonuclear cross-section to use\n"
	    ))
	    .def("set_tear_down_per_call", &I3PropagatorServicePROPOSAL::SetTearDownPerCall)
	;

    class_<PROPOSAL::SimplePropagator, boost::shared_ptr<PROPOSAL::SimplePropagator>,
        boost::noncopyable>(
            "SimplePropagator",
        init<std::string,I3Particle::ParticleType,double,double,double>(
        (arg("medium")="ice",
         arg("type")=I3Particle::MuMinus,
             arg("ecut")=-1.,arg("vcut")=-1,arg("rho")=-1)))
        .def("set_seed", &PROPOSAL::SimplePropagator::SetSeed)
        .def("propagate", &PROPOSAL::SimplePropagator::propagate, (args("p"), arg("distance"), arg("secondaries")=boost::shared_ptr<std::vector<I3Particle> >()))
    ;

}
