#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>

#include "PROPOSAL/PROPOSAL.h"

namespace py = pybind11;
using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Macros
// ------------------------------------------------------------------------- //

#define DEF_PY_PRINT(cls) .def("__str__", &py_print<cls>)

#define COMPONENT_DEF(module, cls)                                                                                     \
    py::class_<Components::cls, Components::Component, std::shared_ptr<Components::cls> >(module, #cls)                \
        .def(py::init<double>(), py::arg("atom_in_molecule") = 1.0);

#define MEDIUM_DEF(module, cls)                                                                                        \
    py::class_<cls, Medium, std::shared_ptr<cls> >(module, #cls)                                                       \
        .def(py::init<double>(), py::arg("density_correction") = 1.0);

#define PARTICLE_DEF(module, cls)                                                                                      \
    py::class_<cls##Def, ParticleDef, std::unique_ptr<cls##Def, py::nodelete> >(module, #cls "Def")                    \
        .def_static("get", &cls##Def::Get, py::return_value_policy::reference);

#define BREMS_DEF(module, cls)                                                                                         \
    py::class_<Brems##cls, std::shared_ptr<Brems##cls>, Bremsstrahlung>(module, #cls)                                  \
        .def(py::init<const ParticleDef&, const Medium&, const EnergyCutSettings&, double, bool>(),                    \
             py::arg("particle_def"),                                                                                  \
             py::arg("medium"),                                                                                        \
             py::arg("energy_cuts"),                                                                                   \
             py::arg("multiplier"),                                                                                    \
             py::arg("lpm_effect"));

#define PHOTO_REAL_DEF(module, cls, parent)                                                                            \
    py::class_<Photo##cls, std::shared_ptr<Photo##cls>, Photo##parent>(module, #cls)                                   \
        .def(py::init<const ParticleDef&, const Medium&, const EnergyCutSettings&, double, bool>(),                    \
             py::arg("particle_def"),                                                                                  \
             py::arg("medium"),                                                                                        \
             py::arg("energy_cuts"),                                                                                   \
             py::arg("multiplier"),                                                                                    \
             py::arg("add_pertubative"));

#define PHOTO_Q2_DEF(module, cls)                                                                                      \
    py::class_<Photo##cls, std::shared_ptr<Photo##cls>, PhotoQ2Integral>(module, #cls)                                 \
        .def(py::init<const ParticleDef&, const Medium&, const EnergyCutSettings&, double, const ShadowEffect&>(),     \
             py::arg("particle_def"),                                                                                  \
             py::arg("medium"),                                                                                        \
             py::arg("energy_cuts"),                                                                                   \
             py::arg("multiplier"),                                                                                    \
             py::arg("shadow_effect"));

#define PHOTO_Q2_INTERPOL_DEF(module, cls)                                                                             \
    py::class_<PhotoQ2Interpolant<Photo##cls>, std::shared_ptr<PhotoQ2Interpolant<Photo##cls> >, Photo##cls>(          \
        module, #cls "Interpolant")                                                                                    \
        .def(py::init<const ParticleDef&,                                                                              \
                      const Medium&,                                                                                   \
                      const EnergyCutSettings&,                                                                        \
                      double,                                                                                          \
                      const ShadowEffect&,                                                                             \
                      InterpolationDef>(),                                                                             \
             py::arg("particle_def"),                                                                                  \
             py::arg("medium"),                                                                                        \
             py::arg("energy_cuts"),                                                                                   \
             py::arg("multiplier"),                                                                                    \
             py::arg("shadow_effect"),                                                                                 \
             py::arg("interpolation_def"));

#define EPAIR_DEF(module, cls)                                                                                         \
    py::class_<Epair##cls, std::shared_ptr<Epair##cls>, EpairProductionRhoIntegral>(module, #cls)                      \
        .def(py::init<const ParticleDef&, const Medium&, const EnergyCutSettings&, double, bool>(),                    \
             py::arg("particle_def"),                                                                                  \
             py::arg("medium"),                                                                                        \
             py::arg("energy_cuts"),                                                                                   \
             py::arg("multiplier"),                                                                                    \
             py::arg("lpm_effect"));

#define EPAIR_INTERPOL_DEF(module, cls)                                                                                \
    py::class_<EpairProductionRhoInterpolant<Epair##cls>,                                                              \
               std::shared_ptr<EpairProductionRhoInterpolant<Epair##cls> >,                                            \
               Epair##cls>(module, #cls "Interpolant")                                                                 \
        .def(py::init<const ParticleDef&, const Medium&, const EnergyCutSettings&, double, bool, InterpolationDef>(),  \
             py::arg("particle_def"),                                                                                  \
             py::arg("medium"),                                                                                        \
             py::arg("energy_cuts"),                                                                                   \
             py::arg("multiplier"),                                                                                    \
             py::arg("lpm_effect"),                                                                                    \
             py::arg("interpolation_def"));

#define MUPAIR_DEF(module, cls)                                                                                         \
    py::class_<Mupair##cls, std::shared_ptr<Mupair##cls>, MupairProductionRhoIntegral>(module, #cls)                      \
        .def(py::init<const ParticleDef&, const Medium&, const EnergyCutSettings&, double>(),                    \
             py::arg("particle_def"),                                                                                  \
             py::arg("medium"),                                                                                        \
             py::arg("energy_cuts"),                                                                                   \
             py::arg("multiplier"));

#define MUPAIR_INTERPOL_DEF(module, cls)                                                                                \
    py::class_<MupairProductionRhoInterpolant<Mupair##cls>,                                                              \
               std::shared_ptr<MupairProductionRhoInterpolant<Mupair##cls> >,                                            \
               Mupair##cls>(module, #cls "Interpolant")                                                                 \
        .def(py::init<const ParticleDef&, const Medium&, const EnergyCutSettings&, double, InterpolationDef>(),  \
             py::arg("particle_def"),                                                                                  \
             py::arg("medium"),                                                                                        \
             py::arg("energy_cuts"),                                                                                   \
             py::arg("multiplier"),                                                                                      \
             py::arg("interpolation_def"));

// ------------------------------------------------------------------------- //
// For __str__
// ------------------------------------------------------------------------- //

template <class T>
std::string py_print(const T& t)
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}

// ------------------------------------------------------------------------- //
// Medium
// ------------------------------------------------------------------------- //

void init_components(py::module& m)
{
    py::module m_sub = m.def_submodule("component");

    m_sub.doc() = R"pbdoc(
        You could create a new component or select one of the implemented.
        Components are used to define a medium as part of a sector.

        A existing component can be called for example with:

        >>> hydro = pyPROPOSAL.component.Hydrogen()
        >>> hydro.atomic_number
        1.00794

        There listed components can be created without initialization:

        * Hydrogen
        * Carbon
        * Nitrogen
        * Oxygen
        * Sodium
        * Magnesium
        * Sulfur
        * Argon
        * Potassium
        * Calcium
        * Iron
        * Copper
        * Lead
        * Uranium
        * StandardRock
        * FrejusRock

        Otherwise you have to initalize the component yourself.
    )pbdoc";

    py::class_<Components::Component, std::shared_ptr<Components::Component>>(m_sub, "Component")
        .def(py::init<std::string, double, double, double>(), py::arg("name"), py::arg("charge"), py::arg("atomic_num"), py::arg("atom_in_molecule"), R"pbdoc(
            Creating a new static component.

            Args:
               name (str):  The name of component.
               charge (float):  Charge in units of Coulomb.
               atomic_num (float): Atom number in periodic table.
               atomic_in_molecule (float): Number of atoms in molecule.
        )pbdoc")
        .def("__str__", &py_print<Components::Component>)
        .def_property_readonly(
                "name", 
                &Components::Component::GetName, 
                R"pbdoc(
                    Get name of component.

                    Returns:
                        str: Name of component
                )pbdoc"
                )
        .def_property_readonly(
                "nuclear_charge",
                &Components::Component::GetNucCharge,
                R"pbdoc(
                    Get nuclear charge of component.

                    Returns:
                        float: Nuclear charge of component
                )pbdoc"
        )
        .def_property_readonly(
                "atomic_number",
                &Components::Component::GetAtomicNum,
                R"pbdoc(
                    Get atomic number of component.

                    Returns:
                        float: Atomic number of component
                )pbdoc"
                )
        .def_property_readonly(
                "atoms_in_molecule",
                &Components::Component::GetAtomInMolecule,
                R"pbdoc(
                    Get number of atoms in one molecule.

                    Returns:
                        float: number of atoms in molecule
                )pbdoc"
                )
        .def_property_readonly(
                "log_constant",
                &Components::Component::GetLogConstant,
                R"pbdoc(
                    Explanation still has to be added.
                )pbdoc"
                )
        .def_property_readonly(
                "bprime",
                &Components::Component::GetBPrime,
                R"pbdoc(
                    Explanation still has to be added.
                )pbdoc"
                )
        .def_property_readonly(
                "average_nucleon_weight",
                &Components::Component::GetAverageNucleonWeight,
                R"pbdoc(
                    Explanation still has to be added.
                )pbdoc"
                )
        .def_property_readonly(
                "mn",
                &Components::Component::GetMN,
                R"pbdoc(
                    Explanation still has to be added.
                )pbdoc"
                )
        .def_property_readonly(
                "r0",
                &Components::Component::GetR0,
                R"pbdoc(
                    Explanation still has to be added.
                )pbdoc"
                );

    COMPONENT_DEF(m_sub, Hydrogen)
    COMPONENT_DEF(m_sub, Carbon)
    COMPONENT_DEF(m_sub, Nitrogen)
    COMPONENT_DEF(m_sub, Oxygen)
    COMPONENT_DEF(m_sub, Sodium)
    COMPONENT_DEF(m_sub, Magnesium)
    COMPONENT_DEF(m_sub, Sulfur)
    COMPONENT_DEF(m_sub, Argon)
    COMPONENT_DEF(m_sub, Potassium)
    COMPONENT_DEF(m_sub, Calcium)
    COMPONENT_DEF(m_sub, Iron)
    COMPONENT_DEF(m_sub, Copper)
    COMPONENT_DEF(m_sub, Lead)
    COMPONENT_DEF(m_sub, Uranium)
    COMPONENT_DEF(m_sub, StandardRock)
    COMPONENT_DEF(m_sub, FrejusRock)
}

void init_medium(py::module& m)
{
    py::module m_sub = m.def_submodule("medium");

    py::class_<Medium, std::shared_ptr<Medium>>(m_sub, "Medium")
        .def("__str__", &py_print<Medium>)
        .def(py::init<>())
        .def(py::init<std::string,
                      double,
                      double,
                      double,
                      double,
                      double,
                      double,
                      double,
                      double,
                      double,
                      const std::vector<std::shared_ptr<Components::Component>>&>(),
             py::arg("name"),
             py::arg("rho"),
             py::arg("I"),
             py::arg("C"),
             py::arg("a"),
             py::arg("m"),
             py::arg("X0"),
             py::arg("X1"),
             py::arg("d0"),
             py::arg("massDensity"),
             py::arg("components"))
        .def_property_readonly("sum_charge", &Medium::GetSumCharge)
        .def_property_readonly("ratio_ZA", &Medium::GetZA)
        .def_property_readonly("ionization_potential", &Medium::GetI)
        .def_property_readonly("refraction_index", &Medium::GetR)
        .def_property_readonly("density_correction", &Medium::GetDensityCorrection)
        .def_property_readonly("radiation_length", &Medium::GetRadiationLength)
        .def_property_readonly("mol_density", &Medium::GetMolDensity)
        .def_property_readonly("average_nucleon_weigth", &Medium::GetMM)
        .def_property_readonly("sum_nucleons", &Medium::GetSumNucleons)
        .def_property_readonly("num_components", &Medium::GetNumComponents)
        .def_property_readonly("components", &Medium::GetComponents)
        .def_property_readonly("name", &Medium::GetName);

    MEDIUM_DEF(m_sub, Water)
    MEDIUM_DEF(m_sub, Ice)
    MEDIUM_DEF(m_sub, Salt)
    MEDIUM_DEF(m_sub, CalciumCarbonate)
    MEDIUM_DEF(m_sub, StandardRock)
    MEDIUM_DEF(m_sub, FrejusRock)
    MEDIUM_DEF(m_sub, Iron)
    MEDIUM_DEF(m_sub, Hydrogen)
    MEDIUM_DEF(m_sub, Lead)
    MEDIUM_DEF(m_sub, Copper)
    MEDIUM_DEF(m_sub, Uranium)
    MEDIUM_DEF(m_sub, Air)
    MEDIUM_DEF(m_sub, Paraffin)
    MEDIUM_DEF(m_sub, AntaresWater)
}

void init_particle(py::module& m)
{
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

    py::class_<ParticleDef, std::shared_ptr<ParticleDef> >(m_sub, "ParticleDef")
        .def(py::init<>())
        .def(py::init<std::string,
                      double,
                      double,
                      double,
                      double,
                      const HardComponentTables::VecType&,
                      const DecayTable&>(),
             py::arg("name"),
             py::arg("mass"),
             py::arg("low"),
             py::arg("lifetime"),
             py::arg("charge"),
             py::arg("hard_component"),
             py::arg("decay_table"))
        .def(py::init<const ParticleDef&>())

        .def("__str__", &py_print<ParticleDef>)
        .def("__eq__", &ParticleDef::operator==)
        .def("__ne__", &ParticleDef::operator!=)
        .def_readonly(
            "name", 
            &ParticleDef::name,
            R"pbdoc(
                name ot the particle
            )pbdoc"
        )
        .def_readonly(
            "mass", 
            &ParticleDef::mass,
            R"pbdoc(
                mass of the particle in MeV
            )pbdoc"
        )
        .def_readonly(
            "low", 
            &ParticleDef::low,
            R"pbdoc(
                lower stochastic sampling energy limit.
            )pbdoc"
        )
        .def_readonly(
            "charge",
             &ParticleDef::charge,
            R"pbdoc(
                charge of the particle
            )pbdoc"
        )
        .def_readonly(
            "decay_table",
             &ParticleDef::decay_table,
            R"pbdoc(
                generate a static particle with in ParticleDefBuilder properties
            )pbdoc"
        )
        // .def_readonly("harc_component_table", &ParticleDef::hard_component_table)
        // .add_property("hard_component_table", make_function(&get_hard_component, return_internal_reference<>()))
        // //TODO(mario): shit Fri 2017/10/13
        ;

    py::class_<ParticleDef::Builder, std::shared_ptr<ParticleDef::Builder>>(m_sub, "ParticleDefBuilder", 
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
            )pbdoc"
        )
        .def(
            py::init<>(),
            R"pbdoc(
                Before you can create or modify a particle a definition builder has 
                to be initalized. It collect the properties of the new or changed 
                particle.
            )pbdoc"
        )
        .def(   
            "SetName", 
            &ParticleDef::Builder::SetName,
            R"pbdoc(
                Args:
                    arg1 (str): name of the particle
            )pbdoc"
        )
        .def(
            "SetMass", 
            &ParticleDef::Builder::SetMass,
            R"pbdoc(
                Args:
                    arg1 (float): mass of the particle
            )pbdoc"
        )
        .def(
            "SetLow", 
            &ParticleDef::Builder::SetLow,
            R"pbdoc(
                Args:
                    arg1 (float): lower stochastic sampling energy limit.
            )pbdoc"
        )
        .def(
            "SetLifetime", 
            &ParticleDef::Builder::SetLifetime,
            R"pbdoc(
                Args:
                    arg1 (float): lifetime of the particle
            )pbdoc"
        )
        .def(
            "SetCharge", 
            &ParticleDef::Builder::SetCharge,
            R"pbdoc(
                Args:
                    arg1 (float): charge of the particle in units of coulomb
            )pbdoc"
        )
        .def(
            "SetDecayTable", 
            &ParticleDef::Builder::SetDecayTable,
            R"pbdoc(
                Args:
                    arg1 (???): ???
            )pbdoc"
        )
        .def(
            "SetParticleDef", 
            &ParticleDef::Builder::SetParticleDef,
            R"pbdoc(
                Args:
                    arg1 (ParticleDef): a pre defined particle which values should 
                        be take over.
            )pbdoc"
        )
        .def(
            "build",
             &ParticleDef::Builder::build,
            R"pbdoc(
                Return: 
                    ParticleDef: generate a static particle with in ParticleDefBuilder properties
            )pbdoc"
        );

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

    py::enum_<DynamicData::Type>(m_sub, "Data")
        .value("None", DynamicData::None)
        .value("Particle", DynamicData::Particle)
        .value("Brems", DynamicData::Brems)
        .value("DeltaE", DynamicData::DeltaE)
        .value("Epair", DynamicData::Epair)
        .value("NuclInt", DynamicData::NuclInt)
        .value("MuPair", DynamicData::MuPair)
        .value("Hadrons", DynamicData::Hadrons)
        .value("ContinuousEnergyLoss", DynamicData::ContinuousEnergyLoss)
        .value("WeakInt", DynamicData::WeakInt);

    py::class_<DynamicData, std::shared_ptr<DynamicData>>(
            m_sub, 
            "DynamicData",
            R"pbdoc(
                Interaction will be stored in form of Dynamic Data. 
                It is used as an array with all important values for 
                secondary particles. Secondary particles are not propagated.
            )pbdoc"
        )
        .def(py::init<DynamicData::Type>())
        .def(py::init<const DynamicData&>())
        .def("__str__", &py_print<DynamicData>)
        .def_property_readonly(
            "id", 
            &DynamicData::GetTypeId,
            R"pbdoc(
                Type of Interaction. Interaction id of a particle can be convertet 
                in an str with:
                
                >>> if p.id == pyPROPOSAL.particle.Data.Particle:
                >>>     print(p.particle_def.name)
                >>> else:
                >>>     print(str(p.id).split(".")[1])
            )pbdoc"
        )
        .def_property(
            "position", 
            &DynamicData::GetPosition, 
            &DynamicData::SetPosition,
            R"pbdoc(
                Place of Interaction.
            )pbdoc"
        )
        .def_property(
            "direction", 
            &DynamicData::GetDirection, 
            &DynamicData::SetDirection,
            R"pbdoc(
                Direction of particle after interaction.
            )pbdoc"
        )
        .def_property(
            "energy", 
            &DynamicData::GetEnergy, 
            &DynamicData::SetEnergy,
            R"pbdoc(
                Energy of secondary particle.
            )pbdoc"
        )
        .def_property(
            "parent_particle_energy", 
            &DynamicData::GetParentParticleEnergy, 
            &DynamicData::SetParentParticleEnergy,
            R"pbdoc(
                Energy of primary particle after interaction.
            )pbdoc"
        )
        .def_property(
            "time", 
            &DynamicData::GetTime, 
            &DynamicData::SetTime,
            R"pbdoc(
                Time since beginning of propagation the primary particle.
            )pbdoc"
        )
        .def_property(
            "propagated_distance", 
            &DynamicData::GetPropagatedDistance, 
            &DynamicData::SetPropagatedDistance,
            R"pbdoc(
                Propagated distance of primary particle.
            )pbdoc"
        );

    py::class_<Particle, std::shared_ptr<Particle>, DynamicData>(
            m_sub, 
            "Particle",
            R"pbdoc(
                The particle class is used as a container to store data
                while propagation process. There every information about
                the primary particle will be stored.

                Information about secondary particles will be found in 
                :meth:`particle.DynamicData`
            )pbdoc"
        )
        .def(py::init<>())
        .def(py::init<const ParticleDef&>())
        .def(py::init<const Particle&>())
        .def(
            "inject_state", 
            &Particle::InjectState,
            R"pbdoc(
            )pbdoc"
        )
        .def_property_readonly(
            "particle_def", 
            &Particle::GetParticleDef,
            R"pbdoc(
                Definition of particle actuell in container
            )pbdoc"
        )
        .def_property_readonly(
            "decay_table", 
            &Particle::GetDecayTable,
            R"pbdoc(
                Decay table of actuell particle
            )pbdoc"
        )
        .def_property(
            "momentum", 
            &Particle::GetMomentum, 
            &Particle::SetMomentum,
            R"pbdoc(
                Momentum of primary particle in eV
            )pbdoc"
        )
        .def_property(
            "entry_point", 
            &Particle::GetEntryPoint, 
            &Particle::SetEntryPoint,
            R"pbdoc(
                Entry point in detector in form of a Vector3d.
            )pbdoc"
        )
        .def_property(
            "entry_time", 
            &Particle::GetEntryTime, 
            &Particle::SetEntryTime, 
            R"pbdoc(
                Time primary particle entered the detector.
            )pbdoc"
        )
        .def_property(
            "entry_energy", 
            &Particle::GetEntryEnergy, 
            &Particle::SetEntryEnergy,
            R"pbdoc(
                Energy primary particle entered the detector.
            )pbdoc"
        )
        .def_property(
            "exit_point",
            &Particle::GetExitPoint, 
            &Particle::SetExitPoint,
            R"pbdoc(
                Point particle exit the detector in form of a Vector3d.
            )pbdoc"
        )
        .def_property(
            "exit_time", 
            &Particle::GetExitTime, 
            &Particle::SetExitTime,
            R"pbdoc(
                Time primary particle exit the detector.
            )pbdoc"
        )
        .def_property(
            "exit_energy", 
            &Particle::GetExitEnergy, 
            &Particle::SetExitEnergy,
            R"pbdoc(
                Energy primary particle exit the detector.
            )pbdoc"
        )
        .def_property(
            "closet_approach_point", 
            &Particle::GetClosestApproachPoint, 
            &Particle::SetClosestApproachPoint,
            R"pbdoc(
                In a first order the point where distance between particle 
                and detector center is minimal.
            )pbdoc" 
        )
        .def_property(
            "closet_approach_time", 
            &Particle::GetClosestApproachTime, 
            &Particle::SetClosestApproachTime,
            R"pbdoc(
                In a first order the time where distance between particle 
                and detector center is minimal.
            )pbdoc"
        )
        .def_property(
            "closet_approach_energy", 
            &Particle::GetClosestApproachEnergy, 
            &Particle::SetClosestApproachEnergy,
            R"pbdoc(
                In a first order the energy where distance between particle 
                and detector center is minimal.
            )pbdoc"
        )
        .def_property(
            "e_lost", 
            &Particle::GetElost, 
            &Particle::SetElost,
            R"pbdoc(
                Energy primary particle lost in detector.
                Energy primary particle lost in detector...
            )pbdoc"
        );
}

void init_decay(py::module& m)
{
    py::module m_sub = m.def_submodule("decay");

    py::class_<DecayChannel, std::shared_ptr<DecayChannel>>(m_sub, "DecayChannel")
        .def("__str__", &py_print<DecayChannel>)
        .def("__eq__", &DecayChannel::operator==)
        .def("__ne__", &DecayChannel::operator!=)
        .def("decay", &DecayChannel::Decay, "Decay the given particle")
        .def_static("boost", (void (*)(Particle&, const Vector3D&, double, double)) &DecayChannel::Boost, "Boost the particle along a direction");

    py::class_<LeptonicDecayChannelApprox, std::shared_ptr<LeptonicDecayChannelApprox>, DecayChannel>(m_sub, "LeptonicDecayChannelApprox")
        .def(py::init<const ParticleDef&, const ParticleDef&, const ParticleDef&>());

    py::class_<LeptonicDecayChannel, std::shared_ptr<LeptonicDecayChannel>, DecayChannel>(m_sub, "LeptonicDecayChannel")
        .def(py::init<const ParticleDef&, const ParticleDef&, const ParticleDef&>());

    py::class_<TwoBodyPhaseSpace, std::shared_ptr<TwoBodyPhaseSpace>, DecayChannel>(m_sub, "TwoBodyPhaseSpace")
        .def(py::init<ParticleDef, ParticleDef>());

    py::class_<ManyBodyPhaseSpace, std::shared_ptr<ManyBodyPhaseSpace>, DecayChannel>(m_sub, "ManyBodyPhaseSpace")
        .def(py::init<std::vector<const ParticleDef*>, PROPOSAL::ManyBodyPhaseSpace::MatrixElementFunction>(),
             py::arg("particle_defs"),
             py::arg("matrix_element") = nullptr)
        .def_static("default_evaluate", &ManyBodyPhaseSpace::DefaultEvaluate, "Return the default matrix element (default 1)")
        .def("evaluate", &ManyBodyPhaseSpace::Evaluate, "Return the matrix element (default 1)")
        .def("set_uniform_sampling", &DecayChannel::SetUniformSampling, "Decide to use uniform phase space sampling");

    py::class_<StableChannel, std::shared_ptr<StableChannel>, DecayChannel>(m_sub, "StableChannel")
        .def(py::init<>());

    py::class_<DecayTable, std::shared_ptr<DecayTable> >(m_sub, "DecayTable")
        .def(py::init<>())
        .def(py::init<const DecayTable&>())
        .def("__str__", &py_print<DecayTable>)
        .def("add_channel", &DecayTable::addChannel, "Add an decay channel")
        .def("select_channel", &DecayTable::SelectChannel, "Select an decay channel according to given branching ratios")
        .def("set_stable", &DecayTable::SetStable, "Define decay table for stable particles")
        .def("set_uniform_sampling", &DecayTable::SetUniformSampling, "Set whether to sample many body decays uniform in phase space");
}

void init_geometry(py::module& m)
{
    py::module m_sub = m.def_submodule("geometry");

    m_sub.doc() = R"pbdoc(
        Every sector is defined by a specific medium and a geometry.
        There are three different classes defined to build a mathematical 
        body. All of them are a object of type :meth:`Geometry`.
        Besides the information of the shape an object of the class 
        geometry contains the position relativ to the coordinate origin.

        Based on the position of the geometry, distances of the propagated 
        particle to geometry sizes can be determined.
    )pbdoc";

    py::enum_<GeometryFactory::Enum>(m_sub, "Shape")
        .value("Sphere", GeometryFactory::Sphere)
        .value("Box", GeometryFactory::Box)
        .value("Cylinder", GeometryFactory::Cylinder);

    py::class_<GeometryFactory::Definition, std::shared_ptr<GeometryFactory::Definition>>(m_sub, "GeometryDefinition")
        .def(py::init<>())
        .def_readwrite(
            "shape", 
            &GeometryFactory::Definition::shape,
            R"pbdoc(
                type of shape of the geometry.
            )pbdoc"
        )
        .def_readwrite(
            "position", 
            &GeometryFactory::Definition::position,
            R"pbdoc(
                position relativ to the coordinates origin.
            )pbdoc"
        )
        .def_readwrite(
            "inner_radius", 
            &GeometryFactory::Definition::inner_radius,
            R"pbdoc(
                inner radius of type :meth:`Sphere`
            )pbdoc"
        )
        .def_readwrite(
            "outer_radius", 
            &GeometryFactory::Definition::radius,
            R"pbdoc(
                inner radius of type :meth:`Sphere`
            )pbdoc"
        )
        .def_readwrite(
            "width", 
            &GeometryFactory::Definition::width,
            R"pbdoc(
                width of type :meth:`Box`
            )pbdoc"
        )
        .def_readwrite(
            "height", 
            &GeometryFactory::Definition::height,
            R"pbdoc(
                height of type :meth:`Box`
            )pbdoc"
        )
        .def_readwrite(
            "depth", 
            &GeometryFactory::Definition::depth,
            R"pbdoc(
                depth of type :meth:`Box`
            )pbdoc"
        );

    py::class_<Geometry, std::shared_ptr<Geometry>>(m_sub, "Geometry")
        DEF_PY_PRINT(Geometry)
        .def(
            "is_infront", 
            &Geometry::IsInfront,
            R"pbdoc(
                Check if particle is in fron of the geometry.

                Parameters:
                    arg0 (Vector3D): particle position
                    arg1 (Vector3D): particle direction

                Return:
                    bool: Is particle in front of the geometry?
            )pbdoc"
        )
        .def(
            "is_inside", 
            &Geometry::IsInside,
            R"pbdoc(
                Check if particle is in the geometry.

                Parameters:
                    arg0 (Vector3D): particle position
                    arg1 (Vector3D): particle direction

                Return:
                    bool: Is particle in the geometry?
            )pbdoc"
        )
        .def(
            "is_behind", 
            &Geometry::IsBehind,
            R"pbdoc(
                Check if particle has passed the geometry.

                Parameters:
                    arg0 (Vector3D): particle position
                    arg1 (Vector3D): particle direction

                Return:
                    bool: Is particle in front of the geometry?
            )pbdoc"
        )
        .def(
            "distance_to_border", 
            &Geometry::DistanceToBorder,
            py::arg("position"),
            py::arg("direction"),
            R"pbdoc(
                Calculates in dependence of the particle position and direction 
                the distance to the next border 

                Parameters:
                    position (Vector3D): particle position
                    direction (Vector3D): particle direction

                Return:
                    float: distance to border
            )pbdoc"
        )
        .def(
            "distance_to_closet_approach", 
            &Geometry::DistanceToClosestApproach,
            R"pbdoc(
                Calculates in dependence of the particle position and direction 
                the distance where the particle pass the geometry center with 
                minimal distance.

                Parameters:
                    arg0 (Vector3D): particle position
                    arg1 (Vector3D): particle direction

                Return:
                    float: distance to closest approach
            )pbdoc"
        )
        .def_property_readonly(
            "name", 
            &Geometry::GetName,
            R"pbdoc(
                name of the geometry
            )pbdoc"
        )
        .def_property(
            "position", 
            &Geometry::GetPosition, 
            &Geometry::SetPosition,
            R"pbdoc(
                position of the geometry center.
            )pbdoc"
        )
        .def_property(
            "hierarchy", 
            &Geometry::GetHierarchy, 
            &Geometry::SetHierarchy,
            R"pbdoc(
                hierachy of the geometry. If sectors overlap, the sector 
                with the highest hierachy will be selected.
            )pbdoc"
        );

    py::class_<Sphere, std::shared_ptr<Sphere>, Geometry>(m_sub, "Sphere")
        .def(py::init<>())
        .def(py::init<Vector3D, double, double>())
        .def(py::init<const Sphere&>())
        .def_property(
            "inner_radius", 
            &Sphere::GetInnerRadius, 
            &Sphere::SetInnerRadius,
            R"pbdoc(
                inner radius of the sphere
            )pbdoc"
        )
        .def_property(
            "radius", 
            &Sphere::GetRadius, 
            &Sphere::SetRadius,
            R"pbdoc(
                outer radius of the sphere
            )pbdoc"
        );

    py::class_<Box, std::shared_ptr<Box>, Geometry>(m_sub, "Box")
        .def(py::init<>())
        .def(py::init<Vector3D, double, double, double>())
        .def(py::init<const Box&>())
        .def_property(
            "width", 
            &Box::GetX, 
            &Box::SetX,
            R"pbdoc(
                width of the box (x-axis)
            )pbdoc"
        )
        .def_property(
            "height", 
            &Box::GetY, 
            &Box::SetY,
            R"pbdoc(
                height of the box (y-axis)
            )pbdoc"
        )
        .def_property(
            "depth", 
            &Box::GetZ, 
            &Box::SetZ,
            R"pbdoc(
                depth of the box (z-axis)
            )pbdoc"
            );

    py::class_<Cylinder, std::shared_ptr<Cylinder>, Geometry>(
            m_sub, 
            "Cylinder",
            R"pbdoc(
                A cylinder can be created as a hollow cylinder.
                For this purpose, a corresponding radius must be 
                selected for the bore along the main axis. A 
                cylinder without a bore is equal to a bore radius 
                equal to zero.
            )pbdoc"
        )
        .def(py::init<>())
        .def(py::init<Vector3D, double, double, double>())
        .def(py::init<const Cylinder&>())
        .def_property(
            "inner_radius", 
            &Cylinder::GetInnerRadius, 
            &Cylinder::SetInnerRadius,
            R"pbdoc(
                the inner radius of the bore through the main axis
                of the cylinder
            )pbdoc"
        )
        .def_property(
            "radius", 
            &Cylinder::GetRadius, 
            &Cylinder::SetRadius,
            R"pbdoc(
                radius of outer shell of the cylinder 
            )pbdoc"
        )
        .def_property(
            "height", &Cylinder::GetZ, 
            &Cylinder::SetZ,
            R"pbdoc(
                height of the cylinder
            )pbdoc"
        );
}

void init_parametrization(py::module& m)
{
    py::module m_sub = m.def_submodule("parametrization");

    // py::class_<Parametrization::IntegralLimits, std::shared_ptr<Parametrization::IntegralLimits>>("IntegralLimits")
    py::class_<Parametrization::IntegralLimits, std::shared_ptr<Parametrization::IntegralLimits>>(m_sub, "IntegralLimits")
        .def(py::init<>())
        .def_readwrite("v_max", &Parametrization::IntegralLimits::vMax,
            R"pbdoc(
            Highest physical possible v for the current parametrization.
            )pbdoc")
        .def_readwrite("v_up", &Parametrization::IntegralLimits::vUp,
            R"pbdoc(
            Energy cut set by the user via the cut settings. Can be energy dependent. Used to differentiate between continous and stochastic losses.

            See :meth:`~pyPROPOSAL.EnergyCutSettings` for more information on the energy cut settings.
            )pbdoc")
        .def_readwrite("v_min", &Parametrization::IntegralLimits::vMin,
            R"pbdoc(
            Lowest physical possible v for the current parametrization
            )pbdoc");

    py::class_<Parametrization, std::shared_ptr<Parametrization>>(m_sub, "Parametrization",
            R"pbdoc(
            Parametrization objects provide the theoretical input for physical cross section used in PROPOSAL, whereas :meth:`~pyPROPOSAL.crosssection.CrossSection`
            provides the numerical methods to process the parametrization. 

            For each physical process in PROPOSAL there are several different parametrizations available, so the user can check how the theoretical input influences
            the simulation.
            )pbdoc")
        .def("__str__", &py_print<Parametrization>)
        .def("differential_crosssection", &Parametrization::DifferentialCrossSection,
         py::arg("energy"),
         py::arg("v"),
                R"pbdoc(
            Calculate the value 

            .. math:: 

                \frac{d\sigma}{dv}(E), 

            e.g. the differential crosssection in v for the current parametrization.
            If the parametrization is given in a double-differential-crosssection the crosssection will be integrated over all values but v.

            Args:
                energy (float): energy in MeV
                v (float): fraction of energy that is transfered from the primary particle to secondary particles (or other forms of energy losses) 

            Return:
                differential_crosssection (float): returns the differential crosssection in v

                )pbdoc"  )
        .def("dEdx_integrand", &Parametrization::FunctionToDEdxIntegral,
         py::arg("energy"),
         py::arg("v"),
                R"pbdoc(
            Calculate the value 

            .. math:: 

                v \cdot \frac{d\sigma}{dv}(E), 

            e.g. the differential crosssection in v for the current parametrization multipied by v. 
            If the parametrization is given in a double-differential-crosssection the crosssection will be integrated over all values but v.

            Args:
                energy (float): energy in MeV
                v (float): fraction of energy that is transfered from the primary particle to secondary particles (or other forms of energy losses) 

            Return:
                dEdx_integrand (float): returns the differential crosssection in v multiplied by v

            PROPOSAL uses this function internally, for example to calcuate :math:`\langle\frac{dE}{dx}\rangle`

                )pbdoc"  )
        .def("dE2dx_integrand", &Parametrization::FunctionToDE2dxIntegral,
         py::arg("energy"),
         py::arg("v"),
                R"pbdoc(
            Calculate the value 

            .. math:: 

                v^2 \cdot \frac{d\sigma}{dv}(E), 

            e.g. the differential crosssection in v for the current parametrization multipied by :math:`v^2`. 
            If the parametrization is given in a double-differential-crosssection the crosssection will be integrated over all values but v.

            Args:
                energy (float): energy in MeV
                v (float): fraction of energy that is transfered from the primary particle to secondary particles (or other forms of energy losses) 

            Return:
                dE2dx_integrand (float): returns the differential crosssection in v multiplied by :math:`v^2`

            PROPOSAL uses this function internally, for example in the calculation of the countinous randomization. 

                )pbdoc"  )
        .def("dNdx_integrand", &Parametrization::FunctionToDNdxIntegral,
         py::arg("energy"),
         py::arg("v"),
                R"pbdoc(
            Calculate the value 

            .. math:: 

                \frac{d\sigma}{dv}(E), 

            e.g. the differential crosssection in v for the current parametrization. 
            If the parametrization is given in a double-differential-crosssection the crosssection will be integrated over all values but v.

            Args:
                energy (float): energy in MeV
                v (float): fraction of energy that is transfered from the primary particle to secondary particles (or other forms of energy losses) 

            Return:
                dNdx_integrand (float): returns the differential crosssection in v

            PROPOSAL uses this function internally, for example to compare the probabilities of the different possible interactions. 

                )pbdoc"  )
        .def("integral_limits", &Parametrization::GetIntegralLimits,
         py::arg("energy"),
                R"pbdoc(
            Returns:
                integral_limits (:meth:`~pyPROPOSAL.parametrization.IntegralLimits`): returns the integral limits for the given energy and the
                current parametrization
                )pbdoc"  
         )
        .def_property_readonly("name", &Parametrization::GetName,
                R"pbdoc(

            Get name of current parametrization

                )pbdoc"
        )
        .def_property_readonly("particle_def", &Parametrization::GetParticleDef,
                R"pbdoc(

            Get :meth:`~pyPROPOSAL.particle.ParticleDef` used by the parametrization

                )pbdoc"
        )
        .def_property_readonly("medium", &Parametrization::GetMedium,
                R"pbdoc( 

            Get :meth:`~pyPROPOSAL.medium` used by the parametrization

                )pbdoc"
        )
        .def_property_readonly("energy_cuts", &Parametrization::GetEnergyCuts,
                R"pbdoc( 

            Get :meth:`~pyPROPOSAL.EnergyCutSettings` defined by the user for the parametrization

                )pbdoc"
        )
        .def_property_readonly("multiplier", &Parametrization::GetMultiplier,
                R"pbdoc( 

            Get multiplier used for the parametrization

                )pbdoc"
        )
        .def_property_readonly("hash", &Parametrization::GetHash,
                R"pbdoc( 

            Get internal hash corresponding to the current parametrization

                )pbdoc"
        );

    // --------------------------------------------------------------------- //
    // Bremsstrahlung
    // --------------------------------------------------------------------- //

    py::module m_sub_brems = m_sub.def_submodule("bremsstrahlung");
    py::class_<Bremsstrahlung, std::shared_ptr<Bremsstrahlung>, Parametrization>(m_sub_brems, "Bremsstrahlung",
                R"pbdoc( 

            Virtual class for the Bremsstrahlung parametrizations. They can be initialized by using one of the given parametrizations with the following parameters

            Args:                                                                                                  
                particle_def (:meth:`~pyPROPOSAL.particle.ParticleDef`): includes all static particle information for the parametrization such as mass, charge, etc.
                medium (:meth:`~pyPROPOSAL.medium`): includes all medium information for the parametrization such as densities or nucleon charges                                   
                energy_cuts (:meth:`~pyPROPOSAL.EnergyCutSettings`): energy cut setting for the parametrization
                multiplier (double): Use a multiplicative factor for the differential crosssection. Can be used for testing or other studies                                                            
                lpm_effect (bool): Enable or disable the corrections due to the Ter-Mikaelian and Landau-Pomeranchuk effect.                                                                                            

            The following parametrizations are currently implemented:

            * KelnerKokoulinPetrukhin

            * PetrukhinShestakov

            * CompleteScreening

            * AndreevBezrukovBugaev

            * SandrockSoedingreksoRhode

            Example:
                To create a bremsstrahlung parametrization

                >>> mu = pyPROPOSAL.particle.MuMinusDef.get()
                >>> medium = pyPROPOSAL.medium.StandardRock(1.0)
                >>> cuts = pyPROPOSAL.EnergyCutSettings(-1, -1)
                >>> param = pyPROPOSAL.parametrization.bremsstrahlung.SandrockSoedingreksoRhode(mu, medium, cuts, 1.0, False)
                )pbdoc"
        );

    BREMS_DEF(m_sub_brems, KelnerKokoulinPetrukhin)
    BREMS_DEF(m_sub_brems, PetrukhinShestakov)
    BREMS_DEF(m_sub_brems, CompleteScreening)
    BREMS_DEF(m_sub_brems, AndreevBezrukovBugaev)
    BREMS_DEF(m_sub_brems, SandrockSoedingreksoRhode)

    py::enum_<BremsstrahlungFactory::Enum>(m_sub_brems, "BremsParametrization")
        .value("PetrukhinShestakov", BremsstrahlungFactory::PetrukhinShestakov)
        .value("KelnerKokoulinPetrukhin", BremsstrahlungFactory::KelnerKokoulinPetrukhin)
        .value("CompleteScreening", BremsstrahlungFactory::CompleteScreening)
        .value("AndreevBezrukovBugaev", BremsstrahlungFactory::AndreevBezrukovBugaev)
        .value("SandrockSoedingreksoRhode", BremsstrahlungFactory::SandrockSoedingreksoRhode);

    py::class_<BremsstrahlungFactory::Definition, std::shared_ptr<BremsstrahlungFactory::Definition> >(m_sub_brems, "BremsDefinition")
        .def(py::init<>())
        .def_readwrite("parametrization", &BremsstrahlungFactory::Definition::parametrization)
        .def_readwrite("lpm_effect", &BremsstrahlungFactory::Definition::lpm_effect)
        .def_readwrite("multiplier", &BremsstrahlungFactory::Definition::multiplier);

    // --------------------------------------------------------------------- //
    // Epair
    // --------------------------------------------------------------------- //

    py::module m_sub_epair = m_sub.def_submodule("pairproduction");
    py::class_<EpairProduction, std::shared_ptr<EpairProduction>, Parametrization>(m_sub_epair, "EpairProduction",
                R"pbdoc( 

            Virtual class for the electron pair production parametrizations. They can be initialized by using one of the given parametrizations with the following parameters

            Args:                                                                                                  
                particle_def (:meth:`~pyPROPOSAL.particle.ParticleDef`): includes all static particle information for the parametrization such as mass, charge, etc.
                medium (:meth:`~pyPROPOSAL.medium`): includes all medium information for the parametrization such as densities or nucleon charges                                   
                energy_cuts (:meth:`~pyPROPOSAL.EnergyCutSettings`): energy cut setting for the parametrization
                multiplier (double): Use a multiplicative factor for the differential crosssection. Can be used for testing or other studies                                                            
                lpm_effect (bool): Enable or disable the corrections due to the Ter-Mikaelian and Landau-Pomeranchuk effect.  
                interpolation_def (:meth:`~pyPROPOSAL.InterpolationDef`): Only needed by Interpolant parametrizations. Includes settings for the interpolation                                    
                                                                                          
            Since the differential cross section is given in :math:`\rho` as well, an intergration over this parameter is needed.
            When using the interpolation_def parameter, this integration is saved in interpolation tables (improving the performance of the calculation with neglible decline in accuracy).

            The following parametrizations are currently implemented:

            * KelnerKokoulinPetrukhin

            * SandrockSoedingreksoRhode

            * KelnerKokoulinPetrukhinInterpolant

            * SandrockSoedingreksoRhodeInterpolant

            Example:
                To create a electron pair production parametrization

                >>> mu = pyPROPOSAL.particle.MuMinusDef.get()
                >>> medium = pyPROPOSAL.medium.StandardRock(1.0)
                >>> cuts = pyPROPOSAL.EnergyCutSettings(-1, -1)
                >>> param = pyPROPOSAL.parametrization.pairproduction.SandrockSoedingreksoRhode(mu, medium, cuts, 1.0, False)
                )pbdoc");

    py::class_<EpairProductionRhoIntegral, std::shared_ptr<EpairProductionRhoIntegral>, EpairProduction>(m_sub_epair, "EpairProductionRhoIntegral")
        .def("function_to_integral", &EpairProductionRhoIntegral::FunctionToIntegral);

    EPAIR_DEF(m_sub_epair, KelnerKokoulinPetrukhin)
    EPAIR_DEF(m_sub_epair, SandrockSoedingreksoRhode)

    EPAIR_INTERPOL_DEF(m_sub_epair, KelnerKokoulinPetrukhin)
    EPAIR_INTERPOL_DEF(m_sub_epair, SandrockSoedingreksoRhode)

    py::enum_<EpairProductionFactory::Enum>(m_sub_epair, "EpairParametrization")
        .value("KelnerKokoulinPetrukhin", EpairProductionFactory::KelnerKokoulinPetrukhin)
        .value("SandrockSoedingreksoRhode", EpairProductionFactory::SandrockSoedingreksoRhode);

    py::class_<EpairProductionFactory::Definition, std::shared_ptr<EpairProductionFactory::Definition> >(m_sub_epair, "EpairDefinition")
        .def(py::init<>())
        .def_readwrite("parametrization", &EpairProductionFactory::Definition::parametrization)
        .def_readwrite("lpm_effect", &EpairProductionFactory::Definition::lpm_effect)
        .def_readwrite("multiplier", &EpairProductionFactory::Definition::multiplier);

    // --------------------------------------------------------------------- //
    // Mupair
    // --------------------------------------------------------------------- //

    py::module m_sub_mupair = m_sub.def_submodule("mupairproduction");
    py::class_<MupairProduction, std::shared_ptr<MupairProduction>, Parametrization>(m_sub_mupair, "MupairProduction",
                R"pbdoc( 

            Virtual class for the muon pair production parametrizations. They can be initialized by using one of the given parametrizations with the following parameters

            Args:                                                                                                  
                particle_def (:meth:`~pyPROPOSAL.particle.ParticleDef`): includes all static particle information for the parametrization such as mass, charge, etc.
                medium (:meth:`~pyPROPOSAL.medium`): includes all medium information for the parametrization such as densities or nucleon charges                                   
                energy_cuts (:meth:`~pyPROPOSAL.EnergyCutSettings`): energy cut setting for the parametrization
                multiplier (double): Use a multiplicative factor for the differential crosssection. Can be used for testing or other studies   
                interpolation_def (:meth:`~pyPROPOSAL.InterpolationDef`): Only needed by Interpolant parametrizations. Includes settings for the interpolation                                    
                                                                                          
            Since the differential cross section is given in :math:`\rho` as well, an intergration over this parameter is needed.
            When using the interpolation_def parameter, this integration is saved in interpolation tables (improving the performance of the calculation with neglible decline in accuracy).                                                         

            The following parametrizations are currently implemented:

            * KelnerKokoulinPetrukhin

            * KelnerKokoulinPetrukhinInterpolant

            Example:
                To create a muon pair production parametrization

                >>> mu = pyPROPOSAL.particle.MuMinusDef.get()
                >>> medium = pyPROPOSAL.medium.StandardRock(1.0)
                >>> cuts = pyPROPOSAL.EnergyCutSettings(-1, -1)
                >>> param = pyPROPOSAL.parametrization.mupairproduction.KelnerKokoulinPetrukhin(mu, medium, cuts, 1.0)
                )pbdoc")
            .def("Calculaterho", &MupairProduction::Calculaterho,
            py::arg("energy"),
            py::arg("v"),
            py::arg("rnd1"),
            py::arg("rnd2"),
                R"pbdoc(
            Muon pairproduction creates a muonpair from the initial particle. The asymmetry between the created muon and antimuon is described by
            the asymmetry parameter rho via 

            .. math:: 

                \rho = \frac{E_+ - E_-}{E_+ + E_-}, 

            with :math:`E_{\pm}` the energy of the created (anti)muon.
            This function samples a value of rho according to the given parameters.

            Args:
                energy (float): energy of initial particle in MeV
                v (float): fraction of energy that is transfered from the primary particle to the created muon pair
                rnd1 (float): random number for sampling the absolute value of rho
                rnd2 (float): random number for sampling the sign of rho 

            Return:
                rho (float): asymmetry parameter :math:`\rho` with :math:`0 \leq \rho \leq 1`

                )pbdoc"  );

    py::class_<MupairProductionRhoIntegral, std::shared_ptr<MupairProductionRhoIntegral>, MupairProduction>(m_sub_mupair, "MupairProductionRhoIntegral")
        .def("function_to_integral", &MupairProductionRhoIntegral::FunctionToIntegral);

    MUPAIR_DEF(m_sub_mupair, KelnerKokoulinPetrukhin)

    MUPAIR_INTERPOL_DEF(m_sub_mupair, KelnerKokoulinPetrukhin)

    py::enum_<MupairProductionFactory::Enum>(m_sub_mupair, "MupairParametrization")
        .value("KelnerKokoulinPetrukhin", MupairProductionFactory::KelnerKokoulinPetrukhin);

    py::class_<MupairProductionFactory::Definition, std::shared_ptr<MupairProductionFactory::Definition> >(m_sub_mupair, "MupairDefinition")
        .def(py::init<>())
        .def_readwrite("parametrization", &MupairProductionFactory::Definition::parametrization)
        .def_readwrite("mupair_enable", &MupairProductionFactory::Definition::mupair_enable)
        .def_readwrite("multiplier", &MupairProductionFactory::Definition::multiplier)
        .def_readwrite("particle_output", &MupairProductionFactory::Definition::particle_output);

    // --------------------------------------------------------------------- //
    // Weak Interaction
    // --------------------------------------------------------------------- //

    py::module m_sub_weak = m_sub.def_submodule("weakinteraction");
    py::class_<WeakInteraction, std::shared_ptr<WeakInteraction>, Parametrization>(m_sub_weak, "WeakInteraction",
                R"pbdoc( 

            Virtual class for the weak interaction parametrizations. They can be initialized by using one of the given parametrizations with the following parameters

            Args:                                                                                                  
                particle_def (:meth:`~pyPROPOSAL.particle.ParticleDef`): includes all static particle information for the parametrization such as mass, charge, etc.
                medium (:meth:`~pyPROPOSAL.medium`): includes all medium information for the parametrization such as densities or nucleon charges                                   
                multiplier (double): Use a multiplicative factor for the differential crosssection. Can be used for testing or other studies                                                            

            The following parametrizations are currently implemented:

            * WeakCooperSarkarMertsch

            Example:
                To create a weak interaction parametrization

                >>> mu = pyPROPOSAL.particle.MuMinusDef.get()
                >>> medium = pyPROPOSAL.medium.StandardRock(1.0)
                >>> param = pyPROPOSAL.parametrization.weakinteraction.WeakCooperSarkarMertsch(mu, medium, 1.0)
                )pbdoc");


    py::class_<WeakCooperSarkarMertsch, std::shared_ptr<WeakCooperSarkarMertsch>, WeakInteraction>(m_sub_weak, "CooperSarkarMertsch")                      \
        .def(py::init<const ParticleDef&, const Medium&, double>(),                    \
             py::arg("particle_def"),                                                                                  \
             py::arg("medium"),                                                                                        \
             py::arg("multiplier"));


    py::enum_<WeakInteractionFactory::Enum>(m_sub_weak, "WeakParametrization")
            .value("CooperSarkarMertsch", WeakInteractionFactory::CooperSarkarMertsch);

    py::class_<WeakInteractionFactory::Definition, std::shared_ptr<WeakInteractionFactory::Definition> >(m_sub_weak, "WeakDefinition")
            .def(py::init<>())
            .def_readwrite("parametrization", &WeakInteractionFactory::Definition::parametrization)
            .def_readwrite("weak_enable", &WeakInteractionFactory::Definition::weak_enable)
            .def_readwrite("multiplier", &WeakInteractionFactory::Definition::multiplier);

    // --------------------------------------------------------------------- //
    // Photo
    // --------------------------------------------------------------------- //

    py::module m_sub_photo = m_sub.def_submodule("photonuclear");
    py::class_<Photonuclear, std::shared_ptr<Photonuclear>, Parametrization>(m_sub_photo, "Photonuclear");

    // Shadow Effect
    py::class_<ShadowEffect, std::shared_ptr<ShadowEffect>>(m_sub_photo, "ShadowEffect",
                R"pbdoc( 

            Virtual class for the parametrizations of the ShadowEffect used in the photonuclear interaction calculations.
            The nucleon shadowing describes the difference in the crosssection between the interaction of a photon with the whole nucleon compared to the interaction of a photon
            with a single nucleon (multplied by the number of nucleons in the atom). In general, the latter cross section in bigger than the real, measured crosssection, therefore
            this effect is called shadowing.
                                                     
            The following parametrizations are currently implemented:

            * ShadowButkevichMikhailov

            * ShadowDuttaRenoSarcevicSeckel

                )pbdoc")
        .def("calculate_shadow_effect", &ShadowEffect::CalculateShadowEffect,
                R"pbdoc( 

            Calculate the shadow effect independently

            Args:                                                                                                  
                Component (:meth:`~pyPROPOSAL.component`): Component to calculate the shadow effect
                x (float): Bjorken x                                
                nu (float): Fraction of energy transfered from the particle via the photon
                                                            
                )pbdoc")
        .def_property_readonly("name", &ShadowEffect::GetName,
                R"pbdoc( 
            Return the name of the current parametrization of the shadow effect                                     
                )pbdoc");

    py::class_<ShadowDuttaRenoSarcevicSeckel, std::shared_ptr<ShadowDuttaRenoSarcevicSeckel>, ShadowEffect>(
        m_sub_photo, "ShadowDuttaRenoSarcevicSeckel")
        .def(py::init<>());

    py::class_<ShadowButkevichMikhailov, std::shared_ptr<ShadowButkevichMikhailov>, ShadowEffect>(
        m_sub_photo, "ShadowButkevichMikhailov")
        .def(py::init<>());

    // Real Photon
    py::class_<RealPhoton, std::shared_ptr<RealPhoton>>(m_sub_photo, "RealPhoton")
        .def("calculate_hard_component", &RealPhoton::CalculateHardComponent)
        .def_property_readonly("name", &RealPhoton::GetName);

    py::class_<SoftComponent, std::shared_ptr<SoftComponent>, RealPhoton>(m_sub_photo, "SoftComponent")
        .def(py::init<>());
    py::class_<HardComponent, std::shared_ptr<HardComponent>>(m_sub_photo,
        "HardComponent")
        .def(py::init<const ParticleDef&>(),py::arg("particle_def"));

    py::class_<PhotoRealPhotonAssumption, std::shared_ptr<PhotoRealPhotonAssumption>, Photonuclear>(m_sub_photo, "PhotoRealPhotonAssumption",
                R"pbdoc( 

            Virtual class for the parametrizations of photonuclear interaction. They can be initialized by using one of the given parametrizations with the following parameters

            Args:                                                                                                  
                particle_def (:meth:`~pyPROPOSAL.particle.ParticleDef`): includes all static particle information for the parametrization such as mass, charge, etc.
                medium (:meth:`~pyPROPOSAL.medium`): includes all medium information for the parametrization such as densities or nucleon charges                                   
                energy_cuts (:meth:`~pyPROPOSAL.EnergyCutSettings`): energy cut setting for the parametrization
                multiplier (double): Use a multiplicative factor for the differential crosssection. Can be used for testing or other studies               
                hard_component (bool): Enabling or disabling the calculation of the hard component
                                                            
            The following parametrizations are currently implemented:

            * Zeus

            * BezrukovBugaev

            * Kokoulin

            * Rhode

            Example:
                To create a photonuclear parametrization

                >>> mu = pyPROPOSAL.particle.MuMinusDef.get()
                >>> medium = pyPROPOSAL.medium.StandardRock(1.0)
                >>> cuts = pyPROPOSAL.EnergyCutSettings(-1, -1)
                >>> param = pyPROPOSAL.parametrization.photonuclear.Rhode(mu, medium, cuts, 1.0, True)
                )pbdoc");
    py::class_<PhotoQ2Integral, std::shared_ptr<PhotoQ2Integral>, Photonuclear>(m_sub_photo, "PhotoQ2Integral",
                R"pbdoc( 

            Virtual class for the parametrizations of photonuclear interaction. They can be initialized by using one of the given parametrizations with the following parameters

            Args:                                                                                                  
                particle_def (:meth:`~pyPROPOSAL.particle.ParticleDef`): includes all static particle information for the parametrization such as mass, charge, etc.
                medium (:meth:`~pyPROPOSAL.medium`): includes all medium information for the parametrization such as densities or nucleon charges                                   
                energy_cuts (:meth:`~pyPROPOSAL.EnergyCutSettings`): energy cut setting for the parametrization
                multiplier (double): Use a multiplicative factor for the differential crosssection. Can be used for testing or other studies               
                ShadowEffect (:meth:`~pyPROPOSAL.parametrization.photonuclear.ShadowEffect`): Parametrization of the ShadowEffect to be used    
                InterpolationDef (:meth:`~pyPROPOSAL.InterpolationDef`): Only needed by Interpolant parametrizations. Includes settings for the interpolation    
 

                                                            
            This virtual class includes all photonuclear interactions where the differential crossection in given in :math:`Q^2`.
            The following parametrizations are currently implemented:

            * AbramowiczLevinLevyMaor91

            * AbramowiczLevinLevyMaor97

            * ButkevichMikhailov

            * RenoSarcevicSu

            * AbramowiczLevinLevyMaor91Interpolant

            * AbramowiczLevinLevyMaor97Interpolant

            * ButkevichMikhailovInterpolant

            * RenoSarcevicSuInterpolant

            The parametrization with "Interpolant" as a suffix creates an interpolation table for the :math:`Q^2` integration, which improved the perfomance.

            Example:
                To create a photonuclear parametrization

                >>> mu = pyPROPOSAL.particle.MuMinusDef.get()
                >>> medium = pyPROPOSAL.medium.StandardRock(1.0)
                >>> cuts = pyPROPOSAL.EnergyCutSettings(-1, -1)
                >>> shadow = pyPROPOSAL.parametrization.photonuclear.ShadowButkevichMikhailov()
                >>> interpol = pyPROPOSAL.InterpolationDef
                >>> param = pyPROPOSAL.parametrization.photonuclear.RenoSarcevicSuInterpolant(mu, medium, cuts, 1.0, shadow, interpol)
                )pbdoc");

    PHOTO_REAL_DEF(m_sub_photo, Zeus, RealPhotonAssumption)
    PHOTO_REAL_DEF(m_sub_photo, BezrukovBugaev, RealPhotonAssumption)
    PHOTO_REAL_DEF(m_sub_photo, Rhode, RealPhotonAssumption)
    PHOTO_REAL_DEF(m_sub_photo, Kokoulin, BezrukovBugaev) // Kokoulin derives from BezrukovBugaev

    PHOTO_Q2_DEF(m_sub_photo, AbramowiczLevinLevyMaor91)
    PHOTO_Q2_DEF(m_sub_photo, AbramowiczLevinLevyMaor97)
    PHOTO_Q2_DEF(m_sub_photo, ButkevichMikhailov)
    PHOTO_Q2_DEF(m_sub_photo, RenoSarcevicSu)

    PHOTO_Q2_INTERPOL_DEF(m_sub_photo, AbramowiczLevinLevyMaor91)
    PHOTO_Q2_INTERPOL_DEF(m_sub_photo, AbramowiczLevinLevyMaor97)
    PHOTO_Q2_INTERPOL_DEF(m_sub_photo, ButkevichMikhailov)
    PHOTO_Q2_INTERPOL_DEF(m_sub_photo, RenoSarcevicSu)

    py::enum_<PhotonuclearFactory::Enum>(m_sub_photo, "PhotoParametrization")
        .value("Zeus", PhotonuclearFactory::Zeus)
        .value("BezrukovBugaev", PhotonuclearFactory::BezrukovBugaev)
        .value("Rhode", PhotonuclearFactory::Rhode)
        .value("Kokoulin", PhotonuclearFactory::Kokoulin)
        .value("AbramowiczLevinLevyMaor91", PhotonuclearFactory::AbramowiczLevinLevyMaor91)
        .value("AbramowiczLevinLevyMaor97", PhotonuclearFactory::AbramowiczLevinLevyMaor97)
        .value("ButkevichMikhailov", PhotonuclearFactory::ButkevichMikhailov)
        .value("RenoSarcevicSu", PhotonuclearFactory::RenoSarcevicSu);

    py::enum_<PhotonuclearFactory::Shadow>(m_sub_photo, "PhotoShadow")
        .value("DuttaRenoSarcevicSeckel", PhotonuclearFactory::ShadowDuttaRenoSarcevicSeckel)
        .value("ButkevichMikhailov", PhotonuclearFactory::ShadowButkevichMikhailov);

    py::class_<PhotonuclearFactory::Definition, std::shared_ptr<PhotonuclearFactory::Definition> >(m_sub_photo, "PhotoDefinition")
        .def(py::init<>())
        .def_readwrite("parametrization", &PhotonuclearFactory::Definition::parametrization)
        .def_readwrite("shadow", &PhotonuclearFactory::Definition::shadow)
        .def_readwrite("hard_component", &PhotonuclearFactory::Definition::hard_component)
        .def_readwrite("multiplier", &PhotonuclearFactory::Definition::multiplier);

    // --------------------------------------------------------------------- //
    // Ionization
    // --------------------------------------------------------------------- //

    py::module m_sub_ioniz = m_sub.def_submodule("ionization");

    py::class_<Ionization, std::shared_ptr<Ionization>, Parametrization>(m_sub_ioniz, "Ionization")
        .def(py::init<const ParticleDef&, const Medium&, const EnergyCutSettings&, double>(),
             py::arg("particle_def"),
             py::arg("medium"),
             py::arg("energy_cuts"),
             py::arg("multiplier"),
                R"pbdoc( 

            Ionization parametrization. It can be initialized with the following parameters

            Args:                                                                                                  
                particle_def (:meth:`~pyPROPOSAL.particle.ParticleDef`): includes all static particle information for the parametrization such as mass, charge, etc.
                medium (:meth:`~pyPROPOSAL.medium`): includes all medium information for the parametrization such as densities or nucleon charges                                   
                energy_cuts (:meth:`~pyPROPOSAL.EnergyCutSettings`): energy cut setting for the parametrization                                                                                           
                multiplier (double): Use a multiplicative factor for the differential crosssection. Can be used for testing or other studies                                                            

            Example:
                To create a ionization parametrization

                >>> mu = pp.particle.MuMinusDef.get()
                >>> medium = pp.medium.StandardRock(1.0)
                >>> cuts = pp.EnergyCutSettings(-1, -1)
                >>> param = pyPROPOSAL.parametrization.ionization.Ionization(mu, medium, cuts, multiplier)
                )pbdoc");

    py::class_<IonizationFactory::Definition, std::shared_ptr<IonizationFactory::Definition> >(m_sub_ioniz, "IonizationDefinition")
        .def(py::init<>())
        .def_readwrite("multiplier", &IonizationFactory::Definition::multiplier);
}

void init_crosssection(py::module& m)
{
    py::module m_sub = m.def_submodule("crosssection");

    py::class_<CrossSection, std::shared_ptr<CrossSection>>(m_sub, "CrossSection",
                R"pbdoc( 

            Virtual class for crosssections. The crosssection class provides all mathematical methods to process the theoretical, differential crosssections that are
            given by the parametrizations. A cross section class can be initialized with the following parameters

            Args:                                                                                                  
                param (:meth:`~pyPROPOSAL.parametrization`): parametrization for the crosssection, including the chosen theoretical model
                interpolation_def (:meth:`~pyPROPOSAL.InterpolationDef`): Only needed by Interpolant parametrizations. Includes settings for the interpolation                                    
 
            The crosssection class can either work with interpolation tables or with exact intergration for every single calculation.
            Since the usage of interpolation tables can improve the speed of the propagation by several orders of magnitude (with neglible decline in accuracy) it is highly recommended
            to use the interpolation methods.
                                                            
            There are specific crosssection classes for every interaction that can be used. 

            * BremsIntegral / BremsInterpolant

            * EpairIntegral / EpairInterpolant

            * IonizIntegral / IonizInterpolant

            * MupairIntegral / MupairInterpolant

            * PhotoIntegral / PhotoInterpolant

            * WeakIntegral / WeakInterpolant

            Example:
                To create a bremsstrahlung CrossSection

                >>> mu = pyPROPOSAL.particle.MuMinusDef.get()
                >>> medium = pyPROPOSAL.medium.StandardRock(1.0)
                >>> cuts = pyPROPOSAL.EnergyCutSettings(-1, -1)
                >>> interpol = pyPROPOSAL.InterpolationDef
                >>> param = pyPROPOSAL.parametrization.bremsstrahlung.SandrockSoedingreksoRhode(mu, medium, cuts, 1.0, False)
                >>> cross = pyPROPOSAL.crosssection.BremsInterpolant(param, interpol)
                >>> cross.calculate_dEdx(1e6) # exmaple usage of the created crosssection class...
                )pbdoc")
        DEF_PY_PRINT(CrossSection)
        .def("calculate_dEdx", &CrossSection::CalculatedEdx,
            py::arg("energy"),
                R"pbdoc( 

            Calculates the continous energy loss :math:`\langle \frac{dE}{dx} \rangle`, which equals to

                .. math:: \frac{N_A}{A} \cdot E \cdot \int_{v_{min}}^{v_{cut}} v \cdot \frac{d\sigma}{dv} dv

            with the particle energy E, the relative energy loss v and the crosssection :math:`\sigma`. The value :math:`v_{cut}` is the energy cut to differentiate between
            continous and stochastic losses in PROPOSAL, see :meth:`~pyPROPOSAL.EnergyCutSettings` for more information on the energy cuts.

            Args:                                                                                                  
                energy (float): energy in MeV
                                                            
                )pbdoc")
        .def("calculate_dE2dx", &CrossSection::CalculatedE2dx,
            py::arg("energy"),
                R"pbdoc( 

            Calculates the value

                .. math:: \frac{N_A}{A} \cdot E^2 \cdot \int_{v_{min}}^{v_{cut}} v^2 \cdot \frac{d\sigma}{dv} dv

            with the particle energy E, the relative energy loss v and the crosssection :math:`\sigma`. The value :math:`v_{cut}` is the energy cut to differentiate between
            continous and stochastic losses in PROPOSAL, see :meth:`~pyPROPOSAL.EnergyCutSettings` for more information on the energy cuts.

            The value is important for the calculation of the ContinuousRandomization (see :meth:`~pyPROPOSAL.ContinuousRandomizer`) 

            Args:                                                                                                  
                energy (float): energy in MeV
                                                            
                )pbdoc")
        .def("calculate_dNdx", (double (CrossSection::*)(double))&CrossSection::CalculatedNdx,
            py::arg("energy"),
                R"pbdoc( 

            Calculates the total cross section

                .. math:: \frac{N_A}{A} \cdot \int_{v_{cut}}^{v_{max}} \frac{d\sigma}{dv} dv

            with the particle energy E, the relative energy loss v and the crosssection :math:`\sigma`. The value v_{cut} is the energy cut to differentiate between
            continous and stochastic losses in PROPOSAL, see :meth:`~pyPROPOSAL.EnergyCutSettings` for more information on the energy cuts.

            Note that this integral only includes the v values about our cut, therefore this values represents only the total crosssection for the stochastic energy losses. 

            Args:                                                                                                  
                energy (float): energy in MeV
                                                            
                )pbdoc")
        .def("calculate_dNdx_rnd", (double (CrossSection::*)(double, double))&CrossSection::CalculatedNdx,
            py::arg("energy"),
            py::arg("rnd"),

                R"pbdoc( 

            Calculates the total cross section

                .. math:: \frac{N_A}{A} \cdot \int_{v_{cut}}^{v_{max}} \frac{d\sigma}{dv} dv

            with the particle energy E, the relative energy loss v and the crosssection :math:`\sigma`. The value v_{cut} is the energy cut to differentiate between
            continous and stochastic losses in PROPOSAL, see :meth:`~pyPROPOSAL.EnergyCutSettings` for more information on the energy cuts.

            Furthermore, for every component in the medium, a stochatic energy loss in saved based on the random number rnd. The values are saved in the crosssection class and can be used
            by other methods.

            Note that this integral only includes the v values about our cut, therefore this values represents only the total crosssection for the stochastic energy losses. 

            Args:                                                                                                  
                energy (float): energy in MeV
                rnd (float): random number between 0 and 1
                                               
                )pbdoc")
        .def("calculate_stochastic_loss", (double (CrossSection::*)(double, double, double))&CrossSection::CalculateStochasticLoss,
            py::arg("energy"),
            py::arg("rnd1"),
            py::arg("rnd2"),
                R"pbdoc( 

            Samples a stochastic energy loss for a particle of the energy E. 

            Args:                                                                                                  
                energy (float): energy in MeV
                rnd1 (float): random number between 0 and 1, samples the energy loss fraction
                rnd2 (float): random number between 0 and 1, sampled the component where the energy loss in occuring

            Returns:
                sampled energy loss for the current particle in MeV

            With inverse transform sampling, using rnd1, the fraction of the energy loss v is determined from the differential crosssection.
            By comparing the total cross sections for every medium, rnd2 is used to determine the component of the current medium for which the stochatic energy loss is calculated.
                     
                )pbdoc")
        .def("calculate_produced_particles", &CrossSection::CalculateProducedParticles,
            py::arg("energy"),
            py::arg("energy_loss"),
            py::arg("rnd1"),
            py::arg("rnd2"),
                R"pbdoc( 

            Available for the muon pairproduction. Samples the two muons that are produces in pairproduction by sampling the asymmetry parameter :math:`\rho` which describes
            the distribution of the energy loss between the two created muons.

            Args:                                                                                                  
                energy (float): primary particle energy in MeV
                energy_loss (float): energy loss of the primary particle in MeV
                rnd1 (float): random number between 0 and 1 to calculate the asymmetry between the two produced particles
                rnd2 (float): random number between 0 and 1 to calculate the sign of the asymmetry between the two produces particles

            Returns:
                list of created particles with corresponding energies
                     
                )pbdoc")
        .def_property_readonly("id", &CrossSection::GetTypeId,
                R"pbdoc( 

            Internal id of the current interaction, see :meth:`~pyPROPOSAL.DynamicData` for all available id's.
                                                            
                )pbdoc")
        .def_property_readonly("parametrization", &CrossSection::GetParametrization,
                R"pbdoc( 

            Pointer to the current parametrization object
                                                            
                )pbdoc");

    py::class_<CrossSectionIntegral, std::shared_ptr<CrossSectionIntegral>, CrossSection>(m_sub, "CrossSectionIntegral");
    py::class_<CrossSectionInterpolant, std::shared_ptr<CrossSectionInterpolant>, CrossSection>(m_sub, "CrossSectionInterpolant");

    py::class_<BremsIntegral, std::shared_ptr<BremsIntegral>, CrossSectionIntegral>(m_sub, "BremsIntegral")
        .def(py::init<const Bremsstrahlung&>(), py::arg("parametrization"));
    py::class_<EpairIntegral, std::shared_ptr<EpairIntegral>, CrossSectionIntegral>(m_sub, "EpairIntegral")
        .def(py::init<const EpairProduction&>(), py::arg("parametrization"));
    py::class_<PhotoIntegral, std::shared_ptr<PhotoIntegral>, CrossSectionIntegral>(m_sub, "PhotoIntegral")
        .def(py::init<const Photonuclear&>(), py::arg("parametrization"));
    py::class_<IonizIntegral, std::shared_ptr<IonizIntegral>, CrossSectionIntegral>(m_sub, "IonizIntegral")
        .def(py::init<const Ionization&>(), py::arg("parametrization"));
    py::class_<MupairIntegral, std::shared_ptr<MupairIntegral>, CrossSectionIntegral>(m_sub, "MupairIntegral")
        .def(py::init<const MupairProduction&>(), py::arg("parametrization"));
    py::class_<WeakIntegral, std::shared_ptr<WeakIntegral>, CrossSectionIntegral>(m_sub, "WeakIntegral")
            .def(py::init<const WeakInteraction&>(), py::arg("parametrization"));

    py::class_<BremsInterpolant, std::shared_ptr<BremsInterpolant>, CrossSectionInterpolant>(
        m_sub, "BremsInterpolant")
        .def(py::init<const Bremsstrahlung&, InterpolationDef>(), py::arg("parametrization"), py::arg("interpolation_def"));
    py::class_<EpairInterpolant, std::shared_ptr<EpairInterpolant>, CrossSectionInterpolant>(
        m_sub, "EpairInterpolant")
        .def(py::init<const EpairProduction&, InterpolationDef>(), py::arg("parametrization"), py::arg("interpolation_def"));
    py::class_<PhotoInterpolant, std::shared_ptr<PhotoInterpolant>, CrossSectionInterpolant>(
        m_sub, "PhotoInterpolant")
        .def(py::init<const Photonuclear&, InterpolationDef>(), py::arg("parametrization"), py::arg("interpolation_def"));
    py::class_<IonizInterpolant, std::shared_ptr<IonizInterpolant>, CrossSectionInterpolant>(
        m_sub, "IonizInterpolant")
        .def(py::init<const Ionization&, InterpolationDef>(), py::arg("parametrization"), py::arg("interpolation_def"));
    py::class_<MupairInterpolant, std::shared_ptr<MupairInterpolant>, CrossSectionInterpolant>(
        m_sub, "MupairInterpolant")
        .def(py::init<const MupairProduction&, InterpolationDef>(), py::arg("parametrization"), py::arg("interpolation_def"));
    py::class_<WeakInterpolant, std::shared_ptr<WeakInterpolant>, CrossSectionInterpolant>(
            m_sub, "WeakInterpolant")
            .def(py::init<const WeakInteraction&, InterpolationDef>(), py::arg("parametrization"), py::arg("interpolation_def"));
}

void init_scattering(py::module& m)
{
    py::module m_sub = m.def_submodule("scattering");

    py::class_<Scattering, std::shared_ptr<Scattering> >(m_sub, "Scattering")
        .def("scatter", &Scattering::Scatter)
        .def_property_readonly("particle", &Scattering::GetParticle);

    py::class_<ScatteringMoliere, std::shared_ptr<ScatteringMoliere>, Scattering>(m_sub, "Moliere")
        .def(py::init<Particle&, const Medium&>());

    py::class_<ScatteringHighlandIntegral, std::shared_ptr<ScatteringHighlandIntegral>, Scattering>(m_sub, "HighlandIntegral")
        .def(py::init<Particle&, Utility&, InterpolationDef>());

    py::class_<ScatteringHighland, std::shared_ptr<ScatteringHighland>, Scattering>(m_sub, "Highland")
        .def(py::init<Particle&, const Medium&>());

    py::class_<ScatteringNoScattering, std::shared_ptr<ScatteringNoScattering>, Scattering>(m_sub, "NoScattering")
        .def(py::init<Particle&, const Medium&>());

    py::enum_<ScatteringFactory::Enum>(m_sub, "ScatteringModel")
        .value("HighlandIntegral", ScatteringFactory::HighlandIntegral)
        .value("Moliere", ScatteringFactory::Moliere)
        .value("Highland", ScatteringFactory::Highland)
        .value("NoScattering", ScatteringFactory::NoScattering);
}

PYBIND11_MODULE(pyPROPOSAL, m)
{
    m.doc() = R"pbdoc(
        .. currentmodule:: pyPROPOSAL
    )pbdoc";

    init_components(m);
    init_medium(m);
    init_particle(m);
    init_decay(m);
    init_geometry(m);
    init_parametrization(m);
    init_crosssection(m);
    init_scattering(m);

    py::class_<Vector3D, std::shared_ptr<Vector3D> >(m, "Vector3D")
        .def(py::init<>())
        .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"), py::arg("z"))
        .def(py::init<const Vector3D&>())
        DEF_PY_PRINT(Vector3D)
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * float())
        .def(float() * py::self)
        .def(py::self * py::self)
        .def(-py::self)

        .def_property_readonly("x", &Vector3D::GetX)
        .def_property_readonly("y", &Vector3D::GetY)
        .def_property_readonly("z", &Vector3D::GetZ)
        .def_property_readonly("radius", &Vector3D::GetRadius)
        .def_property_readonly("phi", &Vector3D::GetPhi)
        .def_property_readonly("theta", &Vector3D::GetTheta)

        .def("set_cartesian_coordinates", &Vector3D::SetCartesianCoordinates)
        .def("set_spherical_coordinates", &Vector3D::SetSphericalCoordinates)
        .def("normalize", &Vector3D::normalise)
        .def("magnitude", &Vector3D::magnitude)
        .def("cartesian_from_spherical", &Vector3D::CalculateCartesianFromSpherical)
        .def("spherical_from_cartesian", &Vector3D::CalculateSphericalCoordinates);

    py::class_<EnergyCutSettings, std::shared_ptr<EnergyCutSettings> >(m, "EnergyCutSettings", 
		R"pbdoc(
			Settings for the lower integration limit.
			Losses below the cut will be handeled continously and the 
			other stochasticaly.

			.. math:: 

				\text{cut} = \begin{cases} e_\text{cut} & E * v_\text{cut} \geq e_\text{cut} \\ v_\text{cut} & \, \text{else} \end{cases}
		)pbdoc")
        .def(
                py::init<>(), 
                R"pbdoc(
                    Initialize some standard settings. 
                    The default e_cut = 500 Mev and v_cut = 0.05 will be set.
                )pbdoc"
            )
        .def(
                py::init<double, double>(), 
                py::arg("ecut"), 
                py::arg("vcut"), R"pbdoc(
                    Set the cut values manualy. 
            
                    Args: 
                        ecut (float): static energy cut. 
                        vcut (float): relativ energy cut.
                )pbdoc")
        .def(py::init<const EnergyCutSettings&>())
        .def("__str__", &py_print<EnergyCutSettings>)
        .def_property(
                "ecut", 
                &EnergyCutSettings::GetEcut, 
                &EnergyCutSettings::SetEcut,
                R"pbdoc(
                    Return set e_cut.

                    Returns:
                        float: e_cut
                )pbdoc"
            )
        .def_property(
                "vcut", 
                &EnergyCutSettings::GetVcut, 
                &EnergyCutSettings::SetVcut,
                R"pbdoc(
                    Return set v_cut.

                    Returns:
                        float: v_cut
                )pbdoc"
            )
        .def(
                "get_cut", 
                &EnergyCutSettings::GetCut, 
                R"pbdoc(
                    Return lower Interpolation/Integration limit.

                    Returns:
                        float: cut
                )pbdoc"
            );

    py::class_<InterpolationDef, std::shared_ptr<InterpolationDef>>(m, "InterpolationDef", 
			R"pbdoc(
				The set standard values have been optimized for performance
				and accuracy. They should not be changed any further 
				without a good reason.

				Example:
					For speed savings it makes sense to specify a path to 
					the tables, to reuse build tables if possible.

					Sometimes it is usefull to look in the tables for a check. 
					To do this binary tables can be diabled.

					>>> interpolDef = pp.InterpolationDef()
					>>> interpolDef.do_binary_tables = False
					>>> interpolDef.path_to_tables = "./custom/table/path"
					>>> interpolDef.path_to_tables_readonly = "./custom/table/path"
			)pbdoc")
        .def(py::init<>())
        .def_readwrite(
			"order_of_interpolation", 
			&InterpolationDef::order_of_interpolation,
			R"pbdoc(
				Order of Interpolation.
			)pbdoc"
		)
        .def_readwrite(
			"path_to_tables", 
			&InterpolationDef::path_to_tables,
			R"pbdoc(
				Path where tables can be written from memory to disk to 
				reuse it if possible.
			)pbdoc"
		)
        .def_readwrite(
			"path_to_tables_readonly", 
			&InterpolationDef::path_to_tables_readonly,
			R"pbdoc(
				Path where tables can be read from disk to avoid to rebuild 
				it.
			)pbdoc"
		)
        .def_readwrite(
			"max_node_energy", 
			&InterpolationDef::max_node_energy,
			R"pbdoc(
				maximum energy that will be interpolated. Energies greater 
				than the value are extrapolated. Default: 1e14 MeV
			)pbdoc"
		)
        .def_readwrite(
			"nodes_cross_section", 
			&InterpolationDef::nodes_cross_section,
			R"pbdoc(
				number of nodes used by evaluation of cross section 
				integrals. Default: xxx
			)pbdoc"
		)
        .def_readwrite(
			"nodes_continous_randomization", 
			&InterpolationDef::nodes_continous_randomization,
			R"pbdoc(
				number of nodes used by evaluation of continous 
				randomization integrals. Default: xxx
			)pbdoc"
		)
        .def_readwrite(
			"nodes_propagate", 
			&InterpolationDef::nodes_propagate,
			R"pbdoc(
				number of nodes used by evaluation of propagation
				integrals. Default: xxx
			)pbdoc"
		)
        .def_readwrite(
			"do_binary_tables", 
			&InterpolationDef::do_binary_tables,
			R"pbdoc(
				Should binary tables be used to store the data. 
				This will increase performance, but are not readable for a 
				crosscheck by human. Default: xxx
			)pbdoc"
		)
        .def_readwrite(
            "just_use_readonly_path", 
            &InterpolationDef::just_use_readonly_path,
            R"pbdoc(
                Just the readonly path to the interpolation tables is used.
                This will stop the program, if the required table is not
                in the readonly path. The (writable) path_to_tables will be
                ignored. Default: xxx
            )pbdoc"
        );

    // --------------------------------------------------------------------- //
    // Utility
    // --------------------------------------------------------------------- //

    py::class_<Utility, std::shared_ptr<Utility> >(m, "Utility")
        .def(py::init<const ParticleDef&, const Medium&, const EnergyCutSettings&, Utility::Definition>(),
             py::arg("partcle_def"),
             py::arg("medium"),
             py::arg("cuts"),
             py::arg("definition"))
        .def(py::init<const ParticleDef&,
                      const Medium&,
                      const EnergyCutSettings&,
                      Utility::Definition,
                      InterpolationDef>(),
             py::arg("partcle_def"),
             py::arg("medium"),
             py::arg("cuts"),
             py::arg("definition"),
             py::arg("interpolation_def"))
        .def_property_readonly("particle_def", &Utility::GetParticleDef)
        .def_property_readonly("medium", &Utility::GetMedium)
        .def_property_readonly("cross_sections", &Utility::GetCrosssections);

    py::class_<Utility::Definition, std::shared_ptr<Utility::Definition> >(m, "UtilityDefinition")
        .def(py::init<>())
        .def_readwrite("brems_def", &Utility::Definition::brems_def)
        .def_readwrite("photo_def", &Utility::Definition::photo_def)
        .def_readwrite("epair_def", &Utility::Definition::epair_def)
        .def_readwrite("ioniz_def", &Utility::Definition::ioniz_def)
        .def_readwrite("mupair_def", &Utility::Definition::mupair_def)
        .def_readwrite("weak_def", &Utility::Definition::weak_def);


    // --------------------------------------------------------------------- //
    // ContinousRandomization
    // --------------------------------------------------------------------- //

    py::class_<ContinuousRandomizer, std::shared_ptr<ContinuousRandomizer> >(m, "ContinuousRandomizer", 
    R"pbdoc(
        If :math:`v_\text{cut}` is large, the spectrum is not continously any 
        more. Particles which has no stochastic loss crossing the medium has 
        all the same energy :math:`E_\text{peak}` after propagating through 
        the medium.  This produce a gap of size of minimal stochastic loss.
        This effect can be reduced using continous randomization.
        
        .. figure:: figures/mu_continuous.png
            :scale: 50%
            :align: center

            Propagate a muon with 100 TeV through 300 m Rock.

        The average energy loss from continous energy losses, is 
        estimated by the energy-integral. 

        .. math::
        
            \int^{E_f}_{E_i} \frac{\sigma(E)}{-f(E)} d\text{E} = -\text{log}(\xi)

        Since probability of a energy loss :math:`e_{lost} < e_{cut}` 
        at distance is finite, it produce fluctuations in average 
        energy losses.

        This small losses can be added in form of a perturbation from 
        average energy loss.

    )pbdoc")
        .def(
                py::init<const Utility&, const InterpolationDef>(),
                py::arg("utility"),
                py::arg("interpolation_def"),
                R"pbdoc(
                    Initalize a continous randomization calculator. 
                    This may take some minutes because for all parametrization 
                    the continous randomization interpolation tables have to be 
                    build.

                    Note:
                        .. math::
                            
                            \langle ( \Delta ( \Delta E ) )^2 \rangle = 
                                \int^{e_\text{cut}}_{e_0} \frac{d\text{E}}{-f(E)} 
                                \left( \int_0^{e_{cut}} e^2 p(e;E) d\text{e} \right)
                    
                    Args:
                        interpolation_def (interpolation_def): specify the number of interpolation points for cont-integral
                        utility (utility): specify the parametrization and energy cuts
                )pbdoc"
            )
        .def(
                "randomize", 
                &ContinuousRandomizer::Randomize,
                py::arg("initial_energy"),
                py::arg("final_energy"),
                py::arg("rand"),
                R"pbdoc(
                    Calculates the stochastical smering of the distribution based on 
                    the second momentum of the parametrizations, the final and intial 
                    energy.
                    
                    Note:
                        A normal distribution with uppere defined variance will be 
                        asumed. The cumulative distibution function has the form:
                        
                        .. math::
                            
                            \text{cdf}(E) = \frac{1}{2} \left(1 + \text{erf} \left(
                                \frac{E}{ \sqrt{ 2 \sigma^2 } } \right) \right)
                        
                        There will be sampled a deviation form mean energy :math:`a`
                        between :math:`\text{cdf}(E_\text{i} - E_\text{f})` and 
                        :math:`\text{cdf}(E_\text{mass} - E_\text{f})` and finaly the 
                        energy updated.

                        .. math::

                            E_\text{f} = \sigma \cdot a + E_\text{f}
                    
                    Args: 
                        initial_energy (float): energy befor stochastical loss 
                        final_energy (float): energy after stochastical loss 
                        rand (float): random number between 0 and 1 

                    Returns:
                        float: randomized final energy
                )pbdoc"
            );

    // --------------------------------------------------------------------- //
    // Sector Definition
    // --------------------------------------------------------------------- //

    // ----[ Location ]-------------------------------------- //

    py::enum_<Sector::ParticleLocation::Enum>(m, "ParticleLocation")
        .value("infront_detector", Sector::ParticleLocation::InfrontDetector)
        .value("inside_detector", Sector::ParticleLocation::InsideDetector)
        .value("behind_detector", Sector::ParticleLocation::BehindDetector);

    py::class_<Sector::Definition, std::shared_ptr<Sector::Definition> >(
            m, 
            "SectorDefinition",
            R"pbdoc(
                Sector definition is a container that collects all important 
                settings for propagation through a sector. There could for 
                example specify the used cross sections or the energy cut 
                settings.

                An example for multiple SectorDefinitions is the standard 
                implementation of propagation. Three sectors are generated,
                which differ in the energy cut settings, to reduce computing 
                power without a measurable loss of accuracy. Infront of the 
                detector energy cut will be like :math:`v_\text{cut} = 0.05`.
                Inside the detector the energy cut :math:`e_\text{cut} = 500 
                \text{MeV}` is chosen so that particles below the detector 
                resolution are statistically sampled. Outside of the detector no
                cuts are made to find a decay point of the particle in first 
                approximation.
            )pbdoc"
            )
        .def(py::init<>())
        .def_readwrite(
                "cut_settings", 
                &Sector::Definition::cut_settings,
                R"pbdoc(
                    Definition of the :meth:`EnergyCutSettings`
                )pbdoc"
        )
        .def_property(
                "medium", 
                &Sector::Definition::GetMedium, 
                &Sector::Definition::SetMedium,
                R"pbdoc(
                    Definition of the :meth:`~pyPROPOSAL.medium.Medium`
                )pbdoc"
        )
        .def_property(
                "geometry", 
                &Sector::Definition::GetGeometry, 
                &Sector::Definition::SetGeometry,
                R"pbdoc(
                    Definiton of the :meth:`~pyPROPOSAL.geometry.Geometry`
                )pbdoc"
        )
        .def_readwrite(
                "do_stochastic_loss_weighting", 
                &Sector::Definition::do_stochastic_loss_weighting,
                R"pbdoc(
                    Boolean value whether the probability of producing a 
                    stochastic loss should be adjusted with a factor, 
                    defaults to False.
                )pbdoc"
        )
        .def_readwrite(
                "stochastic_loss_weighting", 
                &Sector::Definition::stochastic_loss_weighting,
                R"pbdoc(
                    Factor used to scale the probability of producing a 
                    stochastic loss, defaults to 1.0.
                )pbdoc"
        )
        .def_readwrite(
                "stopping_decay", 
                &Sector::Definition::stopping_decay,
                R"pbdoc(
                    
                )pbdoc"
        )
        .def_readwrite(
                "do_continuous_randomization", 
                &Sector::Definition::do_continuous_randomization,
                R"pbdoc(
                    Boolean if continous randomization should be done if 
                    interpolation if is used, defaults to true.
                )pbdoc"
        )
        .def_readwrite(
                "do_continuous_energy_loss_output", 
                &Sector::Definition::do_continuous_energy_loss_output,
                R"pbdoc(
                    
                )pbdoc"
        )
        .def_readwrite(
                "do_exact_time_calculation", 
                &Sector::Definition::do_exact_time_calculation,
                R"pbdoc(
                    Boolean if particle speed could be approach by the speed 
                    of light or should be calculated by solving the track 
                    integral, defaults to false.

                    If the energy is in the order of the rest energy of the 
                    particle, the assumption that the particle moves at the 
                    speed of light becomes increasingly worse. Than it make 
                    sense to calculate the time integral.

                    .. math::
                            
                            t_\text{f} = t_\text{i} + \int_{x_\text{i}}^{x_\text{f}} 
                                \frac{ \text{dx} }{ v(x) }
                )pbdoc"
        )
        .def_readwrite(
                "scattering_model", 
                &Sector::Definition::scattering_model,
                R"pbdoc(
                    Definition of the scattering modell of type :meth:`~pyPROPOSAL.scattering.ScatteringModel`
                    or deactivate scattering.

                    Example:
                        Deactivating scattering can be achieved with:

                        >>> sec = pyPROPOSAL.SectorDefinition()
                        >>> sec.scattering_model = pyPROPOSAL.scattering.ScatteringModel.NoScattering
                )pbdoc"
        )
        .def_readwrite(
                "particle_location", 
                &Sector::Definition::location,
                R"pbdoc(
                    Definition of the relationship of the sectors to each 
                    other of type :meth:`ParticleLocation`.
                )pbdoc"
        )
        .def_readwrite(
                "crosssection_defs", 
                &Sector::Definition::utility_def,
                R"pbdoc(
                    Definition of the crosssection of type :meth:`~pyPROPOSAL.UtilityDefinition`
                )pbdoc"
        );

    // --------------------------------------------------------------------- //
    // Sector
    // --------------------------------------------------------------------- //

    py::class_<Sector, std::shared_ptr<Sector> >(m, "Sector",R"pbdoc(
            A sector is characterized by its homogeneous attitudes. 
            Within a sector there are no boundaries to consider.
		)pbdoc")
        .def(py::init<Particle&, const Sector::Definition&>(), py::arg("particle"), py::arg("sector_definition"))
        .def(py::init<Particle&, const Sector::Definition&, const InterpolationDef&>(),
            py::arg("particle"), py::arg("sector_definition"), py::arg("interpolation_def"))
        .def(
			"propagate", 
			&Sector::Propagate, 
			py::arg("distance"),
			R"pbdoc(
                Args: 
                    distance (float): Distance to propagate in cm.

                Returns:
                    float: if the value is positive, the energy after the propagated distance. If negativ the propagated distance is return with a minus sign.
			)pbdoc"
		)
        .def(
			"CalculateEnergyTillStochastic", 
			&Sector::CalculateEnergyTillStochastic, 
			py::arg("initial_energy"),
			R"pbdoc(
                Samples the energy up to the next stochastic loss. 

                Args:
                    initial_energy (float): Energy in MeV.

                Returns:
                    tuple: (next stochastic energy, decay energy) 
			)pbdoc"
		)
        .def(
			"MakeStochasticLoss", 
			&Sector::MakeStochasticLoss, 
			py::arg("particle_energy"),
			R"pbdoc(
                Samples the stochastic loss.
			
                Args:
                    particle_energy (float): Energy of the propagated 
                        particle.

                Returns:
                    tuple: (energy loss, interaction type)
			)pbdoc"
		)
        .def_property_readonly(
			"particle", 
			&Sector::GetParticle, 
			R"pbdoc(
				Get the internal created particle to modify its properties
				
				Return:
					Particle: propagated Particle
			)pbdoc"
		);

    // --------------------------------------------------------------------- //
    // Randomgenerator
    // --------------------------------------------------------------------- //

    py::class_<RandomGenerator, std::unique_ptr<RandomGenerator, py::nodelete>>(m, "RandomGenerator")
        .def("random_double", &RandomGenerator::RandomDouble)
        .def("set_seed", &RandomGenerator::SetSeed, py::arg("seed") = 0)
        .def("set_random_function", &RandomGenerator::SetRandomNumberGenerator, py::arg("function"))
        .def_static("get", &RandomGenerator::Get, py::return_value_policy::reference);

    // --------------------------------------------------------------------- //
    // Propagator
    // --------------------------------------------------------------------- //

    py::class_<Propagator, std::shared_ptr<Propagator>>(m, "Propagator")
        .def(
                py::init<const ParticleDef&, const std::vector<Sector::Definition>&, const Geometry&>(), 
                py::arg("particle_def"), 
                py::arg("sector_defs"), 
                py::arg("detector")
            )
        .def(
                py::init<const ParticleDef&, const std::vector<Sector::Definition>&, const Geometry&, const InterpolationDef&>(), 
                py::arg("particle_def"), 
                py::arg("sector_defs"), 
                py::arg("detector"), 
                py::arg("interpolation_def"), R"pbdoc(
                    Function Docstring.

                    Args:
                        particle_def (pyPROPOSAL.particle.ParticleDef): sldkfjas lsdkjfa 
                        sector_defs (List[pyPROPOSAL.SectorDefinition]): slkd sldkf ksjdlfa
                        detector (pyPROPOSAL.geometry.Geometry): kfsdljflaskdj
                        interpolation_def (pyPROPOSAL.InterpolationDef):  sdfa  sldkf sdfa
                )pbdoc")
        .def(
                py::init<const ParticleDef&, 
                const std::vector<Sector::Definition>&, const Geometry&>(),
                py::arg("particle_def"), 
                py::arg("sector_defs"), 
                py::arg("detector")
            )
        .def(
                py::init<const ParticleDef&, const std::string&>(), 
                py::arg("particle_def"), 
                py::arg("config_file")
            )
        .def(
                "propagate", 
                &Propagator::Propagate, 
                py::arg("max_distance_cm") = 1e20, 
                py::return_value_policy::reference,
                R"pbdoc(
                    Propagate a particle through sectors and produce stochastic
                    losses, untill propagated distance is reached.

                    Args:
                        max_distance_cm (float): Maximum distance a particle is
                            propagated before it is considered lost.

                    Returns:
                        list(DynamicData): list of stochastic losses parameters
                        list(list): list of stochastic losses parameters
                    
                    Example:
                        Propagate 1000 particle with an inital energy of 100 
                        Tev and save the losses in daughters.

                        >>> for i in range(int(1e3)):
                        >>>   mu.energy = 1e8
                        >>>   mu.propagated_distance = 0
                        >>>   mu.position = pp.Vector3D(0, 0, 0)
                        >>>   mu.direction = pp.Vector3D(0, 0, -1)
                        >>>   daughters = prop.propagate()

                    Further basic condition are set at the next point of 
                    interaction during generation.
                    The aim is to introduce forced stochastic losses
                    so that the particle can be propagated through homogeneous 
                    sectors.
                    A more percise description of propagation through a
                    homoegenous sector can be found in ???.
                    
                    .. figure:: figures/sector.png
                        :height: 200px
                        :align: center

                        If the next interaction point is outside the actuell 
                        sector, the particle is forced to an interaction at the
                        sector boundary. Thus it can be assumed that the 
                        particle is porpagated by a homogeneous medium.

                    The propagated particle will be forced to interact in the 
                    closest point to the dector center. 
                    This will be stored in the closest_approach variable of 
                    the propagated particle.

                    The propagation terminate if the maximal distance is 
                    reached.
                    The energy loss of the particle (e_lost) in the detector 
                    will be calculated and the produced secondary particles 
                    returned.
                )pbdoc"
            )
        .def_property_readonly(
                "particle", 
                &Propagator::GetParticle, 
                R"pbdoc(
                    Get the internal created particle to modify its properties.

                    Returns:
                        ParticleDef: definition of the propagated particle.
                )pbdoc"
            )
        .def_property_readonly(
                "detector", 
                &Propagator::GetDetector, 
                R"pbdoc(
                    Get the detector geometry.

                    Returns:
                        SectorDefinition: definition of the sector.
                )pbdoc"
            );

    // --------------------------------------------------------------------- //
    // PropagatorService
    // --------------------------------------------------------------------- //

    py::class_<PropagatorService, std::shared_ptr<PropagatorService>>(m, "PropagatorService")
        .def(py::init<>())
        .def("propagate", &PropagatorService::Propagate, py::arg("particle"), py::arg("distance") = 1e20)
        .def("register_propagator", &PropagatorService::RegisterPropagator, py::arg("propagator"));
}

#undef COMPONENT_DEF
#undef MEDIUM_DEF
#undef PARTICLE_DEF
#undef BREMS_DEF
#undef PHOTO_REAL_DEF
#undef PHOTO_Q2_DEF
#undef PHOTO_Q2_INTERPOL_DEF
#undef EPAIR_DEF
#undef EPAIR_INTERPOL_DEF
#undef MUPAIR_DEF
#undef MUPAIR_INTERPOL_DEF
