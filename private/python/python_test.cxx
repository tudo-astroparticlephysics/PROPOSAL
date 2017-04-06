#include <boost/python.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <string>
#include "PROPOSAL/PROPOSALParticle.h"
#include "PROPOSAL/Propagator.h"


// BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getName_overloads , PROPOSALParticle::GetName, 0, 1)

BOOST_PYTHON_MODULE(pyPROPOSAL)
{

    using namespace boost::python;

    // --------------------------------------------------------------------- //
    // ParticleType
    // --------------------------------------------------------------------- //


    enum_<PROPOSALParticle::ParticleType>("ParticleType")
            .value("EPlus",                PROPOSALParticle::ParticleType::EPlus)
            .value("EMinus",               PROPOSALParticle::ParticleType::EPlus)
            .value("MuPlus",               PROPOSALParticle::ParticleType::MuPlus)
            .value("MuMinus",              PROPOSALParticle::ParticleType::MuMinus)
            .value("TauPlus",              PROPOSALParticle::ParticleType::TauPlus)
            .value("TauMinus",             PROPOSALParticle::ParticleType::TauMinus)
            .value("NuE",                  PROPOSALParticle::ParticleType::NuE)
            .value("NuEBar",               PROPOSALParticle::ParticleType::NuEBar)
            .value("NuMu",                 PROPOSALParticle::ParticleType::NuMu)
            .value("NuMuBar",              PROPOSALParticle::ParticleType::NuMuBar)
            .value("NuTau",                PROPOSALParticle::ParticleType::NuTau)
            .value("NuTauBar",             PROPOSALParticle::ParticleType::NuTauBar)
            .value("Brems",                PROPOSALParticle::ParticleType::Brems)
            .value("DeltaE",               PROPOSALParticle::ParticleType::DeltaE)
            .value("EPair",                PROPOSALParticle::ParticleType::EPair)
            .value("NuclInt",              PROPOSALParticle::ParticleType::NuclInt)
            .value("MuPair",               PROPOSALParticle::ParticleType::MuPair)
            .value("Hadrons",              PROPOSALParticle::ParticleType::Hadrons)
            .value("ContinuousEnergyLoss", PROPOSALParticle::ParticleType::ContinuousEnergyLoss)
            .value("Monopole",             PROPOSALParticle::ParticleType::Monopole)
            .value("STauPlus",             PROPOSALParticle::ParticleType::STauPlus)
            .value("STauMinus",            PROPOSALParticle::ParticleType::STauMinus)
            ;

    // --------------------------------------------------------------------- //
    // Particle
    // --------------------------------------------------------------------- //

    std::string (PROPOSALParticle::*getName1)() const = &PROPOSALParticle::GetName;

    class_<PROPOSALParticle, boost::shared_ptr<PROPOSALParticle>>("Particle",
                                                                  init<PROPOSALParticle::ParticleType>(
                                                                  (arg("type")=PROPOSALParticle::ParticleType::MuMinus)))

            .def(self_ns::str(self_ns::self))

            .add_property("energy", &PROPOSALParticle::GetEnergy, &PROPOSALParticle::SetEnergy)
            .add_property("propagated_distance", &PROPOSALParticle::GetPropagatedDistance, &PROPOSALParticle::SetPropagatedDistance)
            .add_property("X", &PROPOSALParticle::GetX, &PROPOSALParticle::SetX)
            .add_property("Y", &PROPOSALParticle::GetY, &PROPOSALParticle::SetY)
            .add_property("Z", &PROPOSALParticle::GetY, &PROPOSALParticle::SetZ)
            .add_property("T", &PROPOSALParticle::GetT, &PROPOSALParticle::SetT)
            .add_property("theta", &PROPOSALParticle::GetTheta, &PROPOSALParticle::SetTheta)
            .add_property("phi", &PROPOSALParticle::GetPhi, &PROPOSALParticle::SetPhi)
            .add_property("momentum", &PROPOSALParticle::GetMomentum, &PROPOSALParticle::SetMomentum)
            .add_property("mass", &PROPOSALParticle::GetMass, &PROPOSALParticle::SetMass)
            .add_property("lifetime", &PROPOSALParticle::GetLifetime, &PROPOSALParticle::SetLifetime)
            .add_property("charge", &PROPOSALParticle::GetCharge, &PROPOSALParticle::SetCharge)
            .add_property("name", getName1)
            .add_property("low", &PROPOSALParticle::GetLow, &PROPOSALParticle::SetLow)
            .add_property("type", &PROPOSALParticle::GetType, &PROPOSALParticle::SetType)
            .add_property("parent_particle_id", &PROPOSALParticle::GetParentParticleId, &PROPOSALParticle::SetParentParticleId)
            .add_property("parent_particle_energy", &PROPOSALParticle::GetParentParticleEnergy, &PROPOSALParticle::SetParentParticleEnergy)
            .add_property("particle_id", &PROPOSALParticle::GetParticleId, &PROPOSALParticle::SetParticleId)

            .add_property("Xi", &PROPOSALParticle::GetXi, &PROPOSALParticle::SetXi)
            .add_property("Yi", &PROPOSALParticle::GetYi, &PROPOSALParticle::SetYi)
            .add_property("Zi", &PROPOSALParticle::GetZi, &PROPOSALParticle::SetZi)
            .add_property("Ti", &PROPOSALParticle::GetTi, &PROPOSALParticle::SetTi)
            .add_property("Ei", &PROPOSALParticle::GetEi, &PROPOSALParticle::SetEi)

            .add_property("Xf", &PROPOSALParticle::GetXf, &PROPOSALParticle::SetXf)
            .add_property("Yf", &PROPOSALParticle::GetYf, &PROPOSALParticle::SetYf)
            .add_property("Zf", &PROPOSALParticle::GetZf, &PROPOSALParticle::SetZf)
            .add_property("Tf", &PROPOSALParticle::GetTf, &PROPOSALParticle::SetTf)
            .add_property("Ef", &PROPOSALParticle::GetEf, &PROPOSALParticle::SetEf)

            .add_property("Xc", &PROPOSALParticle::GetXc, &PROPOSALParticle::SetXc)
            .add_property("Yc", &PROPOSALParticle::GetYc, &PROPOSALParticle::SetYc)
            .add_property("Zc", &PROPOSALParticle::GetZc, &PROPOSALParticle::SetZc)
            .add_property("Tc", &PROPOSALParticle::GetTc, &PROPOSALParticle::SetTc)
            .add_property("Ec", &PROPOSALParticle::GetEc, &PROPOSALParticle::SetEc)

            .add_property("energy_lost", &PROPOSALParticle::GetElost, &PROPOSALParticle::SetElost)
        ;

    // --------------------------------------------------------------------- //
    // Propagator
    // --------------------------------------------------------------------- //

    class_<std::vector<PROPOSALParticle*>>("Secondarys")
        .def(vector_indexing_suite<std::vector<PROPOSALParticle*>>())
        ;

    class_<Propagator, boost::shared_ptr<Propagator>>("Propagator",
                                                      init<std::string, PROPOSALParticle*, bool>(
                                                      (arg("config"),
                                                       arg("particle"),
                                                       arg("applyoptions")=true)))

            .def("propagate", &Propagator::propagate, (arg("max_distance_cm") = 1e20))
            .def("apply_options", &Propagator::ApplyOptions)
            .def("reset_particle", &Propagator::ResetParticle)

            // .add_property("particle", &Propagator::SetParticle)
            .add_property("particle", make_function(&Propagator::GetParticle, return_value_policy<reference_existing_object>()), &Propagator::SetParticle)
            .add_property("seed",&Propagator::GetSeed ,&Propagator::SetSeed)
            .add_property("brems",&Propagator::GetBrems ,&Propagator::SetBrems)
            .add_property("photo",&Propagator::GetPhoto ,&Propagator::SetPhoto)
            .add_property("path_to_tables",&Propagator::GetPath_to_tables ,&Propagator::SetPath_to_tables)
            .add_property("stopping_decay",&Propagator::GetStopping_decay ,&Propagator::SetStopping_decay)
        ;
}
