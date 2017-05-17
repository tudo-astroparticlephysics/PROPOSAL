
#include "PROPOSAL-icetray/I3PropagatorServicePROPOSAL.h"
#include "PROPOSAL-icetray/SimplePropagator.h"
#include "icetray/load_project.h"

I3_PYTHON_MODULE(PROPOSAL)
{
    load_project("PROPOSAL", false);

    using namespace PROPOSAL;
    using namespace boost::python;
    import("icecube.dataclasses");
    import("icecube.sim_services");


    enum_<ParametrizationType::Enum>("CrossSectionParametrization")
        .value("BremsKelnerKokoulinPetrukhin"                  , ParametrizationType::BremsKelnerKokoulinPetrukhin)
        .value("BremsPetrukhinShestakov"                       , ParametrizationType::BremsPetrukhinShestakov)
        .value("BremsAndreevBezrrukovBugaev"                   , ParametrizationType::BremsAndreevBezrukovBugaev)
        .value("BremsCompleteScreeningCase"                    , ParametrizationType::BremsCompleteScreeningCase)
        .value("PhotoKokoulinShadowBezrukovSoft"               , ParametrizationType::PhotoKokoulinShadowBezrukovSoft)
        .value("PhotoKokoulinShadowBezrukovHard"               , ParametrizationType::PhotoKokoulinShadowBezrukovHard)
        .value("PhotoRhodeShadowBezrukovSoft"                  , ParametrizationType::PhotoRhodeShadowBezrukovSoft)
        .value("PhotoRhodeShadowBezrukovHard"                  , ParametrizationType::PhotoRhodeShadowBezrukovHard)
        .value("PhotoBezrukovBugaevShadowBezrukovSoft"         , ParametrizationType::PhotoBezrukovBugaevShadowBezrukovSoft)
        .value("PhotoBezrukovBugaevShadowBezrukovHard"         , ParametrizationType::PhotoBezrukovBugaevShadowBezrukovHard)
        .value("PhotoZeusShadowBezrukovSoft"                   , ParametrizationType::PhotoZeusShadowBezrukovSoft)
        .value("PhotoZeusShadowBezrukovHard"                   , ParametrizationType::PhotoZeusShadowBezrukovHard)
        .value("PhotoAbramowiczLevinLevyMaor91ShadowDutta"     , ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowDutta)
        .value("PhotoAbramowiczLevinLevyMaor91ShadowButkevich" , ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowButkevich)
        .value("PhotoAbramowiczLevinLevyMaor97ShadowDutta"     , ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowDutta)
        .value("PhotoAbramowiczLevinLevyMaor97ShadowButkevich" , ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich)
        .value("PhotoButkevichMikhailovShadowDutta"            , ParametrizationType::PhotoButkevichMikhailovShadowDutta)
        .value("PhotoButkevichMikhailovShadowButkevich"        , ParametrizationType::PhotoButkevichMikhailovShadowButkevich)
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
                 arg("photonuclearParametrization")=ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich
                ),
                 ":param mediadef: Path the the media definition file\n"
                 ":param tabledir: Path to a directory in which to store interpolation constants for cross-section integrals\n"
                 ":param cylinderRadius: Radius of the target volume in meters\n"
                 ":param cylinderHeight: Full height of the target volume in meters\n"
                 ":param type: Type of particle to propagate\n"
                 ":param particleMass: Mass of the propagated particle in GeV. This is only used if type is something exotic.\n"
                 ":param bremsstrahlungParametrization: Parametrization of the bremsstrahlung cross-section to use\n"
                 ":param photonuclearParametrization: Parametrization of photonuclear cross-section to use\n"
            )
        )
        .def("set_tear_down_per_call", &I3PropagatorServicePROPOSAL::SetTearDownPerCall)
    ;

    class_<PROPOSAL::SimplePropagator, boost::shared_ptr<PROPOSAL::SimplePropagator>,
        boost::noncopyable>(
            "SimplePropagator",
            init<std::string,I3Particle::ParticleType,double,double,double>(
                (arg("medium")="ice",
                 arg("type")=I3Particle::MuMinus,
                 arg("ecut")=-1.,
                 arg("vcut")=-1,
                 arg("rho")=-1
                )
            )
        )
        .def("set_seed", &PROPOSAL::SimplePropagator::SetSeed)
        .def("propagate", &PROPOSAL::SimplePropagator::propagate, 
             (args("p"), arg("distance"), arg("secondaries")=boost::shared_ptr<std::vector<I3Particle> >())
            )
    ;

}
