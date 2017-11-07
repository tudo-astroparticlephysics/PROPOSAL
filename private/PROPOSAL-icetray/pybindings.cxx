
// #include "PROPOSAL-icetray/I3PropagatorServicePROPOSAL.h"
#include "PROPOSAL-icetray/SimplePropagator.h"
#include "icetray/load_project.h"

I3_PYTHON_MODULE(PROPOSAL)
{
    load_project("PROPOSAL", false);

    using namespace PROPOSAL;
    using namespace boost::python;
    import("icecube.dataclasses");
    import("icecube.sim_services");


    // enum_<ParametrizationType::Enum>("CrossSectionParametrization")
    //     .value("BremsKelnerKokoulinPetrukhin"                  , ParametrizationType::BremsKelnerKokoulinPetrukhin)
    //     .value("BremsPetrukhinShestakov"                       , ParametrizationType::BremsPetrukhinShestakov)
    //     .value("BremsAndreevBezrrukovBugaev"                   , ParametrizationType::BremsAndreevBezrukovBugaev)
    //     .value("BremsCompleteScreeningCase"                    , ParametrizationType::BremsCompleteScreeningCase)
    //     .value("PhotoKokoulinShadowBezrukovSoft"               , ParametrizationType::PhotoKokoulinShadowBezrukovSoft)
    //     .value("PhotoKokoulinShadowBezrukovHard"               , ParametrizationType::PhotoKokoulinShadowBezrukovHard)
    //     .value("PhotoRhodeShadowBezrukovSoft"                  , ParametrizationType::PhotoRhodeShadowBezrukovSoft)
    //     .value("PhotoRhodeShadowBezrukovHard"                  , ParametrizationType::PhotoRhodeShadowBezrukovHard)
    //     .value("PhotoBezrukovBugaevShadowBezrukovSoft"         , ParametrizationType::PhotoBezrukovBugaevShadowBezrukovSoft)
    //     .value("PhotoBezrukovBugaevShadowBezrukovHard"         , ParametrizationType::PhotoBezrukovBugaevShadowBezrukovHard)
    //     .value("PhotoZeusShadowBezrukovSoft"                   , ParametrizationType::PhotoZeusShadowBezrukovSoft)
    //     .value("PhotoZeusShadowBezrukovHard"                   , ParametrizationType::PhotoZeusShadowBezrukovHard)
    //     .value("PhotoAbramowiczLevinLevyMaor91ShadowDutta"     , ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowDutta)
    //     .value("PhotoAbramowiczLevinLevyMaor91ShadowButkevich" , ParametrizationType::PhotoAbramowiczLevinLevyMaor91ShadowButkevich)
    //     .value("PhotoAbramowiczLevinLevyMaor97ShadowDutta"     , ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowDutta)
    //     .value("PhotoAbramowiczLevinLevyMaor97ShadowButkevich" , ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich)
    //     .value("PhotoButkevichMikhailovShadowDutta"            , ParametrizationType::PhotoButkevichMikhailovShadowDutta)
    //     .value("PhotoButkevichMikhailovShadowButkevich"        , ParametrizationType::PhotoButkevichMikhailovShadowButkevich)
    // ;
    //
    // enum_<MediumType::Enum>("Medium")
    //     .value("Water"         , MediumType::Water)
    //     .value("Ice"           , MediumType::Ice)
    //     .value("Hydrogen"      , MediumType::Hydrogen)
    //     .value("Iron"          , MediumType::Iron)
    //     .value("Copper"        , MediumType::Copper)
    //     .value("Lead"          , MediumType::Lead)
    //     .value("Uranium"       , MediumType::Uranium)
    //     .value("Air"           , MediumType::Air)
    //     .value("AntaresWater"  , MediumType::AntaresWater)
    //     .value("StandardRock"  , MediumType::StandardRock)
    //     .value("FrejusRock"    , MediumType::FrejusRock)
    //     .value("Salt"          , MediumType::Salt)
    //     .value("MineralOil"    , MediumType::MineralOil)
    // ;


    // class_<I3PropagatorServicePROPOSAL, boost::shared_ptr<I3PropagatorServicePROPOSAL>,
    //     bases<I3PropagatorService>, boost::noncopyable>(
    //         "I3PropagatorServicePROPOSAL",
    //         init<I3Particle&, std::string>(
    //             (arg("particle"),
    //              arg("config_file")=I3PropagatorServicePROPOSAL::GetDefaultConfigFileDir(),
    //             ),
    //              ":param particle: Parametrization of the bremsstrahlung cross-section to use\n"
    //              ":param config_file: Path to the config file\n"
    //         )
    //     )
    //     // .def("set_tear_down_per_call", &I3PropagatorServicePROPOSAL::SetTearDownPerCall)
    // ;
    //
    // class_<PROPOSAL::SimplePropagator, boost::shared_ptr<PROPOSAL::SimplePropagator>,
    //     boost::noncopyable>(
    //         "SimplePropagator",
    //         init<I3Particle::ParticleType, MediumType::Enum, double, double, double>(
    //             (arg("type")=I3Particle::MuMinus,
    //              arg("medium")="ice",
    //              arg("ecut")=-1.,
    //              arg("vcut")=-1,
    //              arg("rho")=-1
    //             )
    //         )
    //     )
    //     .def("set_seed", &PROPOSAL::SimplePropagator::SetSeed)
    //     .def("propagate", &PROPOSAL::SimplePropagator::propagate,
    //          (args("p"), arg("distance"), arg("secondaries")=boost::shared_ptr<std::vector<I3Particle> >())
    //         )
    // ;

}
