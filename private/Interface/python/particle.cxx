#include <string>

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/Secondaries.h"
#include "pyBindings.h"

#define PARTICLE_DEF(module, cls)                                                    \
    py::class_<cls##Def, ParticleDef, std::shared_ptr<cls##Def>>(module, #cls "Def") \
        .def(py::init<>());                                                          \

namespace py = pybind11;
using namespace PROPOSAL;

void init_particle(py::module& m) {
    py::module m_sub = m.def_submodule("particle");

    m_sub.doc() = R"pbdoc(
        For each propagation a defined particle is needed.
        You have the possibility to define one by your own or select one of
        the listed below. Particle data are based on source ??.

        +----------+----------+----------+----------+----------+----------+
        |MuMinus   |MuPlus    |EMinus    |EPlus     |TauMinus  |TauPlus   |
        +----------+----------+----------+----------+----------+----------+
        |StauMinus |StauPlus  |Pi0       |PiMinus   |PiPlus    |K0        |
        +----------+----------+----------+----------+----------+----------+
        |KMinus    |KPlus     |PMinus    |PPlus     |NuE       |NuEBar    |
        +----------+----------+----------+----------+----------+----------+
        |NuMu      |NuMuBar   |NuTau     |NuTauBar  |Monopole  |Gamma     |
        +----------+----------+----------+----------+----------+----------+
        |SMPMinus  |SMPPlus   |          |          |          |          |
        +----------+----------+----------+----------+----------+----------+

        Particle are static objects, so that the properties can not be
        changed once it is initalized.  Own particle can be created if you
        initalize a complete new with the :meth:`ParticleDef` class or
        modifie a copy of a existing one with the :meth:`ParticleDefBuilder`.

        A predefined particle can be initialize for example with

        >>> muon = proposal.particle.MuMinusDef.get()

        The :meth:`Particle` class is a container for partilce related data.
        There can be for example initial energy set and losses read out.
        The :meth:`particle.Data` class is a container for interaction a
        particle does.
    )pbdoc";

    py::class_<ParticleDef, std::shared_ptr<ParticleDef>>(m_sub, "ParticleDef")
        .def(py::init<>())
        .def(py::init<std::string, double, double, double, double,
                      const HardComponentTables::VecType&, const DecayTable&,
                      const int, const int>(),
             py::arg("name"), py::arg("mass"), py::arg("low"),
             py::arg("lifetime"), py::arg("charge"), py::arg("hard_component"),
             py::arg("decay_table"), py::arg("particle_type"), py::arg("weak_partner"))
        .def(py::init<const ParticleDef&>())

        .def("__str__", &py_print<ParticleDef>)
        .def("__eq__", &ParticleDef::operator==)
        .def("__ne__", &ParticleDef::operator!=)
        .def_readonly("name", &ParticleDef::name,
                      R"pbdoc(
                name ot the particle
            )pbdoc")
        .def_readonly("mass", &ParticleDef::mass,
                      R"pbdoc(
                mass of the particle in MeV
            )pbdoc")
        .def_readonly("low", &ParticleDef::low,
                      R"pbdoc(
                lower stochastic sampling energy limit.
            )pbdoc")
        .def_readonly("lifetime", &ParticleDef::lifetime,
                      R"pbdoc(
                average lifetime of the particle in seconds.
            )pbdoc")
        .def_readonly("charge", &ParticleDef::charge,
                      R"pbdoc(
                charge of the particle
            )pbdoc")
        .def_readonly("decay_table", &ParticleDef::decay_table,
                      R"pbdoc(
                generate a static particle with in ParticleDefBuilder properties
            )pbdoc")
        .def_readonly("particle_type", &ParticleDef::particle_type,
                      R"pbdoc(
                particle type of the particle
            )pbdoc")
        .def_readonly("weak_partner", &ParticleDef::weak_partner,
                      R"pbdoc(
                particle type of the weak partner particle
            )pbdoc")
        // .def_readonly("hard_component_table",
        // &ParticleDef::hard_component_table)
        // .add_property("hard_component_table",
        // make_function(&get_hard_component, return_internal_reference<>()))
        // //TODO(mario): shit Fri 2017/10/13
        ;

    py::class_<ParticleDef::Builder, std::shared_ptr<ParticleDef::Builder>>(
        m_sub, "ParticleDefBuilder",
        R"pbdoc(
                With the ParticleDefBuilder it is possible to define new particles
                or change porperties of existing one. First you have to initalize
                a ParticleDefBuilder. The new properties of the particle can be
                collected before a new static particle is build.

                Example:
                    A usecase for the ParticleDefBuilder might be if your
                    particle will be propagated stochastically until a
                    certain lower limit which is higher than the particle
                    mass is reached.

                    >>> builder = pp.particle.ParticleDefBuilder()
                    >>> builder.SetParticleDef(pp.particle.MuMinusDef.get())
                    >>> builder.SetLow(1e6)  # 1 Tev
                    >>> mu_def = mu_def_builder.build()

                    Therby muons which energy is lower than one TeV will be
                    handled continously.
            )pbdoc")
        .def(py::init<>(),
             R"pbdoc(
                Before you can create or modify a particle a definition builder has
                to be initalized. It collect the properties of the new or changed
                particle.
            )pbdoc")
        .def("SetName", &ParticleDef::Builder::SetName,
             R"pbdoc(
                Args:
                    arg1 (str): name of the particle
            )pbdoc")
        .def("SetMass", &ParticleDef::Builder::SetMass,
             R"pbdoc(
                Args:
                    arg1 (float): mass of the particle
            )pbdoc")
        .def("SetLow", &ParticleDef::Builder::SetLow,
             R"pbdoc(
                Args:
                    arg1 (float): lower stochastic sampling energy limit.
            )pbdoc")
        .def("SetLifetime", &ParticleDef::Builder::SetLifetime,
             R"pbdoc(
                Args:
                    arg1 (float): lifetime of the particle
            )pbdoc")
        .def("SetCharge", &ParticleDef::Builder::SetCharge,
             R"pbdoc(
                Args:
                    arg1 (float): charge of the particle in units of coulomb
            )pbdoc")
        .def("SetDecayTable", &ParticleDef::Builder::SetDecayTable,
             R"pbdoc(
                Args:
                    arg1 (???): ???
            )pbdoc")
        .def("SetParticleType", &ParticleDef::Builder::SetParticleType,
             R"pbdoc(
                Args:
                    arg1 (int): particle type of the particle
            )pbdoc")
        .def("SetWeakPartner", &ParticleDef::Builder::SetWeakPartner,
             R"pbdoc(
                Args:
                    arg1 (ParticleType): particle type of the weak partner particle
            )pbdoc")
        .def("SetParticleDef", &ParticleDef::Builder::SetParticleDef,
             R"pbdoc(
                Args:
                    arg1 (ParticleDef): a pre defined particle which values should
                        be take over.
            )pbdoc")
        .def("build", &ParticleDef::Builder::build,
             R"pbdoc(
                Return:
                    ParticleDef: generate a static particle with in ParticleDefBuilder properties
            )pbdoc");

    PARTICLE_DEF(m_sub, MuMinus)
    PARTICLE_DEF(m_sub, MuPlus)
    PARTICLE_DEF(m_sub, EMinus)
    PARTICLE_DEF(m_sub, EPlus)
    PARTICLE_DEF(m_sub, TauMinus)
    PARTICLE_DEF(m_sub, TauPlus)
    PARTICLE_DEF(m_sub, StauMinus)
    PARTICLE_DEF(m_sub, StauPlus)
    PARTICLE_DEF(m_sub, Pi0)
    PARTICLE_DEF(m_sub, PiMinus)
    PARTICLE_DEF(m_sub, PiPlus)
    PARTICLE_DEF(m_sub, K0)
    PARTICLE_DEF(m_sub, KMinus)
    PARTICLE_DEF(m_sub, KPlus)
    PARTICLE_DEF(m_sub, PMinus)
    PARTICLE_DEF(m_sub, PPlus)
    PARTICLE_DEF(m_sub, NuE)
    PARTICLE_DEF(m_sub, NuEBar)
    PARTICLE_DEF(m_sub, NuMu)
    PARTICLE_DEF(m_sub, NuMuBar)
    PARTICLE_DEF(m_sub, NuTau)
    PARTICLE_DEF(m_sub, NuTauBar)
    PARTICLE_DEF(m_sub, Monopole)
    PARTICLE_DEF(m_sub, Gamma)
    PARTICLE_DEF(m_sub, SMPMinus)
    PARTICLE_DEF(m_sub, SMPPlus)

    py::class_<DynamicData, std::shared_ptr<DynamicData>>(m_sub, "DynamicData",
                                                          R"pbdoc(
                Interaction will be stored in form of Dynamic Data.
                It is used as an array with all important values for
                secondary particles. Secondary particles are not propagated.
            )pbdoc")
        .def(py::init<const int&>())
        .def(py::init<const DynamicData&>())
        .def("__str__", &py_print<DynamicData>)
        .def_property_readonly("type", &DynamicData::GetType,
                               R"pbdoc(
                Type of Dynamic Data. Can be an Interaction Type
                or a Particle Type.
            )pbdoc")
        .def_property_readonly("name", &DynamicData::GetName,
                               R"pbdoc(
                Name of Particle or interaction.
            )pbdoc")
        .def_property("position", &DynamicData::GetPosition,
                      &DynamicData::SetPosition,
                      R"pbdoc(
                Place of Interaction.
            )pbdoc")
        .def_property("direction", &DynamicData::GetDirection,
                      &DynamicData::SetDirection,
                      R"pbdoc(
                Direction of particle after interaction.
            )pbdoc")
        .def_property("energy", &DynamicData::GetEnergy,
                      &DynamicData::SetEnergy,
                      R"pbdoc(
                Energy of secondary particle.
            )pbdoc")
        .def_property("momentum", &DynamicData::GetMomentum,
                      &DynamicData::SetMomentum,
                      R"pbdoc(
                Momentum of primary particle in MeV
            )pbdoc")
        .def_property("parent_particle_energy",
                      &DynamicData::GetParentParticleEnergy,
                      &DynamicData::SetParentParticleEnergy,
                      R"pbdoc(
                Energy of primary particle after interaction.
            )pbdoc")
        .def_property("time", &DynamicData::GetTime, &DynamicData::SetTime,
                      R"pbdoc(
                Time since beginning of propagation the primary particle.
            )pbdoc")
        .def_property("propagated_distance",
                      &DynamicData::GetPropagatedDistance,
                      &DynamicData::SetPropagatedDistance,
                      R"pbdoc(
                Propagated distance of primary particle.
            )pbdoc");

    py::class_<Secondaries, std::shared_ptr<Secondaries>>(m_sub, "Secondaries",
            R"pbdoc(List of secondaries)pbdoc")
        .def("Query", overload_cast_<const int&>()(&Secondaries::Query, py::const_), py::arg("Interaction"))
        .def("Query", overload_cast_<const std::string&>()(&Secondaries::Query, py::const_), py::arg("Interaction"))
        .def("decay", &Secondaries::DoDecay)
        .def_property_readonly("particles", &Secondaries::GetSecondaries)
        .def_property_readonly("number_of_particles", &Secondaries::GetNumberOfParticles)
        .def_property_readonly("position", &Secondaries::GetPosition)
        .def_property_readonly("direction", &Secondaries::GetDirection)
        .def_property_readonly("parent_particle_energy", &Secondaries::GetParentParticleEnergy)
        .def_property_readonly("energy", &Secondaries::GetEnergy)
        .def_property_readonly("time", &Secondaries::GetTime)
        .def_property_readonly("propagated_distance", &Secondaries::GetPropagatedDistance)
        .def_property_readonly("entry_point", &Secondaries::GetEntryPoint)
        .def_property_readonly("exit_point", &Secondaries::GetExitPoint)
        .def_property_readonly("closest_approach_point", &Secondaries::GetClosestApproachPoint);

    py::enum_<InteractionType>(m_sub, "Interaction_Type")
        .value("Particle", InteractionType::Particle)
        .value("Brems", InteractionType::Brems)
        .value("DeltaE", InteractionType::DeltaE)
        .value("Epair", InteractionType::Epair)
        .value("NuclInt", InteractionType::NuclInt)
        .value("MuPair", InteractionType::MuPair)
        .value("Hadrons", InteractionType::Hadrons)
        .value("ContinuousEnergyLoss", InteractionType::ContinuousEnergyLoss)
        .value("Compton", InteractionType::Compton)
        .value("WeakInt", InteractionType::WeakInt);

   py::enum_<ParticleType>(m_sub, "Particle_Type")
        .value("None", ParticleType::None)
        .value("EMinus", ParticleType::EMinus)
        .value("EPlus", ParticleType::EPlus)
        .value("NuE", ParticleType::NuE)
        .value("NuEBar", ParticleType::NuEBar)
        .value("MuMinus", ParticleType::MuMinus)
        .value("NuMu", ParticleType::NuMu)
        .value("NuMuBar", ParticleType::NuMuBar)
        .value("MuPlus", ParticleType::MuPlus)
        .value("TauMinus", ParticleType::TauMinus)
        .value("TauPlus", ParticleType::TauPlus)
        .value("NuTau", ParticleType::NuTau)
        .value("NuTauBar", ParticleType::NuTauBar)
        .value("Gamma", ParticleType::Gamma)
        .value("Pi0", ParticleType::Pi0)
        .value("PiPlus", ParticleType::PiPlus)
        .value("PiMinus", ParticleType::PiMinus)
        .value("K0", ParticleType::K0)
        .value("KPlus", ParticleType::KPlus)
        .value("KMinus", ParticleType::KMinus)
        .value("STauMinus", ParticleType::STauMinus)
        .value("STauPlus", ParticleType::STauPlus)
        .value("PPlus", ParticleType::PPlus)
        .value("PMinus", ParticleType::PMinus)
        .value("Monopole", ParticleType::Monopole);
}

#undef PARTICLE_DEF
