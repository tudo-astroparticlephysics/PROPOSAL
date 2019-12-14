#include <string>

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/Secondaries.h"
#include "pyBindings.h"

#define PARTICLE_DEF(module, cls)                                           \
    py::class_<cls##Def, ParticleDef,                                       \
               std::unique_ptr<cls##Def, py::nodelete>>(module, #cls "Def") \
        .def_static("get", &cls##Def::Get,                                  \
                    py::return_value_policy::reference);

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

        >>> muon = pyPROPOSAL.particle.MuMinusDef.get()

        The :meth:`Particle` class is a container for partilce related data.
        There can be for example initial energy set and losses read out.
        The :meth:`particle.Data` class is a container for interaction a
        particle does.
    )pbdoc";

    py::class_<ParticleDef, std::shared_ptr<ParticleDef>>(m_sub, "ParticleDef")
        .def(py::init<>())
        .def(py::init<std::string, double, double, double, double,
                      const HardComponentTables::VecType&, const DecayTable&, const int>(),
             py::arg("name"), py::arg("mass"), py::arg("low"),
             py::arg("lifetime"), py::arg("charge"), py::arg("hard_component"),
             py::arg("decay_table"), py::arg("particle_type"))
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
        .def_readonly("charge", &ParticleDef::charge,
                      R"pbdoc(
                charge of the particle
            )pbdoc")
        .def_readonly("decay_table", &ParticleDef::decay_table,
                      R"pbdoc(
                generate a static particle with in ParticleDefBuilder properties
            )pbdoc")
        // .def_readonly("harc_component_table",
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

    py::class_<InteractionType, std::shared_ptr<InteractionType>>(m_sub, "interation_type_id")
        .def_readonly_static("Particle", &InteractionType::Particle)
        .def_readonly_static("Brems", &InteractionType::Brems)
        .def_readonly_static("DeltaE", &InteractionType::DeltaE)
        .def_readonly_static("Epair", &InteractionType::Epair)
        .def_readonly_static("NuclInt", &InteractionType::NuclInt)
        .def_readonly_static("MuPair", &InteractionType::MuPair)
        .def_readonly_static("Hadrons", &InteractionType::Hadrons)
        .def_readonly_static("ContinuousEnergyLoss", &InteractionType::ContinuousEnergyLoss)
        .def_readonly_static("Compton", &InteractionType::Compton)
        .def_readonly_static("WeakInt", &InteractionType::WeakInt);

    py::class_<ParticleType, std::shared_ptr<ParticleType>>(m_sub, "particle_type_id")
        .def_readonly_static("None", &ParticleType::None)
        .def_readonly_static("EMinus", &ParticleType::EMinus)
        .def_readonly_static("EPlus", &ParticleType::EPlus)
        .def_readonly_static("NuE", &ParticleType::NuE)
        .def_readonly_static("NuEBar", &ParticleType::NuEBar)
        .def_readonly_static("MuMinus", &ParticleType::MuMinus)
        .def_readonly_static("NuMu", &ParticleType::NuMu)
        .def_readonly_static("NuMuBar", &ParticleType::NuMuBar)
        .def_readonly_static("MuPlus", &ParticleType::MuPlus)
        .def_readonly_static("TauMinus", &ParticleType::TauMinus)
        .def_readonly_static("TauPlus", &ParticleType::TauPlus)
        .def_readonly_static("NuTau", &ParticleType::NuTau)
        .def_readonly_static("NuTauBar", &ParticleType::NuTauBar)
        .def_readonly_static("Gamma", &ParticleType::Gamma)
        .def_readonly_static("Pi0", &ParticleType::Pi0)
        .def_readonly_static("PiPlus", &ParticleType::PiPlus)
        .def_readonly_static("PiMinus", &ParticleType::PiMinus)
        .def_readonly_static("K0", &ParticleType::K0)
        .def_readonly_static("KPlus", &ParticleType::KPlus)
        .def_readonly_static("KMinus", &ParticleType::KMinus)
        .def_readonly_static("STauMinus", &ParticleType::STauMinus)
        .def_readonly_static("STauPlus", &ParticleType::STauPlus)
        .def_readonly_static("PPlus", &ParticleType::PPlus)
        .def_readonly_static("PMinus", &ParticleType::PMinus)
        .def_readonly_static("Monopole", &ParticleType::Monopole);

    py::class_<DynamicData::Type, std::shared_ptr<DynamicData::Type>>(m_sub, "dyndata_type_id")
        .def_readonly_static("Particle", &DynamicData::Type::Particle)
        .def_readonly_static("Brems", &DynamicData::Type::Brems)
        .def_readonly_static("DeltaE", &DynamicData::Type::DeltaE)
        .def_readonly_static("Epair", &DynamicData::Type::Epair)
        .def_readonly_static("NuclInt", &DynamicData::Type::NuclInt)
        .def_readonly_static("MuPair", &DynamicData::Type::MuPair)
        .def_readonly_static("Hadrons", &DynamicData::Type::Hadrons)
        .def_readonly_static("ContinuousEnergyLoss", &DynamicData::Type::ContinuousEnergyLoss)
        .def_readonly_static("Compton", &DynamicData::Type::Compton)
        .def_readonly_static("WeakInt", &DynamicData::Type::WeakInt)
        .def_readonly_static("EMinus", &DynamicData::Type::EMinus)
        .def_readonly_static("EPlus", &DynamicData::Type::EPlus)
        .def_readonly_static("NuE", &DynamicData::Type::NuE)
        .def_readonly_static("NuEBar", &DynamicData::Type::NuEBar)
        .def_readonly_static("MuMinus", &DynamicData::Type::MuMinus)
        .def_readonly_static("NuMu", &DynamicData::Type::NuMu)
        .def_readonly_static("NuMuBar", &DynamicData::Type::NuMuBar)
        .def_readonly_static("MuPlus", &DynamicData::Type::MuPlus)
        .def_readonly_static("TauMinus", &DynamicData::Type::TauMinus)
        .def_readonly_static("TauPlus", &DynamicData::Type::TauPlus)
        .def_readonly_static("NuTau", &DynamicData::Type::NuTau)
        .def_readonly_static("NuTauBar", &DynamicData::Type::NuTauBar)
        .def_readonly_static("Gamma", &DynamicData::Type::Gamma)
        .def_readonly_static("Pi0", &DynamicData::Type::Pi0)
        .def_readonly_static("PiPlus", &DynamicData::Type::PiPlus)
        .def_readonly_static("PiMinus", &DynamicData::Type::PiMinus)
        .def_readonly_static("K0", &DynamicData::Type::K0)
        .def_readonly_static("KPlus", &DynamicData::Type::KPlus)
        .def_readonly_static("KMinus", &DynamicData::Type::KMinus)
        .def_readonly_static("STauMinus", &DynamicData::Type::STauMinus)
        .def_readonly_static("STauPlus", &DynamicData::Type::STauPlus)
        .def_readonly_static("PPlus", &DynamicData::Type::PPlus)
        .def_readonly_static("PMinus", &DynamicData::Type::PMinus)
        .def_readonly_static("Monopole", &DynamicData::Type::Monopole);

    py::class_<DynamicData, std::shared_ptr<DynamicData>>(m_sub, "DynamicData",
                                                          R"pbdoc(
                Interaction will be stored in form of Dynamic Data.
                It is used as an array with all important values for
                secondary particles. Secondary particles are not propagated.
            )pbdoc")
        .def(py::init<const int&>())
        .def(py::init<const DynamicData&>())
        .def("__str__", &py_print<DynamicData>)
        .def_property_readonly("id", &DynamicData::GetTypeId,
                               R"pbdoc(
                Type of Interaction. Interaction id of a particle can be convertet
                in an str with:

                >>> if p.id == pyPROPOSAL.particle.Data.Particle:
                >>>     print(p.particle_def.name)
                >>> else:
                >>>     print(str(p.id).split(".")[1])
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

    py::class_<Particle, std::shared_ptr<Particle>, DynamicData>(m_sub,
                                                                 "Particle",
                                                                 R"pbdoc(
                The particle class is used as a container to store data
                while propagation process. There every information about
                the primary particle will be stored.

                Information about secondary particles will be found in
                :meth:`particle.DynamicData`
            )pbdoc")
        .def(py::init<>())
        .def(py::init<const ParticleDef&>())
        .def(py::init<const Particle&>())
        .def("inject_state", &Particle::InjectState,
             R"pbdoc(
            )pbdoc")
        .def_property_readonly("particle_def", &Particle::GetParticleDef,
                               R"pbdoc(
                Definition of particle actuell in container
            )pbdoc")
        .def_property_readonly("decay_table", &Particle::GetDecayTable,
                               R"pbdoc(
                Decay table of actuell particle
            )pbdoc")
        .def_property("momentum", &Particle::GetMomentum,
                      &Particle::SetMomentum,
                      R"pbdoc(
                Momentum of primary particle in eV
            )pbdoc")
        .def_property("entry_point", &Particle::GetEntryPoint,
                      &Particle::SetEntryPoint,
                      R"pbdoc(
                Entry point in detector in form of a Vector3d.
            )pbdoc")
        .def_property("entry_time", &Particle::GetEntryTime,
                      &Particle::SetEntryTime,
                      R"pbdoc(
                Time primary particle entered the detector.
            )pbdoc")
        .def_property("entry_energy", &Particle::GetEntryEnergy,
                      &Particle::SetEntryEnergy,
                      R"pbdoc(
                Energy primary particle entered the detector.
            )pbdoc")
        .def_property("exit_point", &Particle::GetExitPoint,
                      &Particle::SetExitPoint,
                      R"pbdoc(
                Point particle exit the detector in form of a Vector3d.
            )pbdoc")
        .def_property("exit_time", &Particle::GetExitTime,
                      &Particle::SetExitTime,
                      R"pbdoc(
                Time primary particle exit the detector.
            )pbdoc")
        .def_property("exit_energy", &Particle::GetExitEnergy,
                      &Particle::SetExitEnergy,
                      R"pbdoc(
                Energy primary particle exit the detector.
            )pbdoc")
        .def_property("closet_approach_point",
                      &Particle::GetClosestApproachPoint,
                      &Particle::SetClosestApproachPoint,
                      R"pbdoc(
                In a first order the point where distance between particle
                and detector center is minimal.
            )pbdoc")
        .def_property("closet_approach_time", &Particle::GetClosestApproachTime,
                      &Particle::SetClosestApproachTime,
                      R"pbdoc(
                In a first order the time where distance between particle
                and detector center is minimal.
            )pbdoc")
        .def_property("closet_approach_energy",
                      &Particle::GetClosestApproachEnergy,
                      &Particle::SetClosestApproachEnergy,
                      R"pbdoc(
                In a first order the energy where distance between particle
                and detector center is minimal.
            )pbdoc")
        .def_property("e_lost", &Particle::GetElost, &Particle::SetElost,
                      R"pbdoc(
                Energy primary particle lost in detector.
                Energy primary particle lost in detector...
            )pbdoc");

    py::class_<Secondaries, std::shared_ptr<Secondaries>>(m_sub, "Secondaries",
            R"pbdoc(List of secondaries)pbdoc")
        .def("Query", overload_cast_<const int&>()(&Secondaries::Query, py::const_), py::arg("Interaction"))
        .def("Query", overload_cast_<const std::string&>()(&Secondaries::Query, py::const_), py::arg("Interaction"))
        .def_property_readonly("particles", &Secondaries::GetSecondaries)
        .def_property_readonly("number_of_particles", &Secondaries::GetNumberOfParticles)
        .def_property_readonly("position", &Secondaries::GetDirection)
        .def_property_readonly("direction", &Secondaries::GetPosition)
        .def_property_readonly("parent_particle_energy", &Secondaries::GetParentParticleEnergy)
        .def_property_readonly("energy", &Secondaries::GetEnergy)
        .def_property_readonly("time", &Secondaries::GetTime)
        .def_property_readonly("propagated_distance", &Secondaries::GetPropagatedDistance);

}

#undef PARTICLE_DEF
