
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
        .def(py::init<const ParticleDef&, const Medium&, const EnergyCutSettings&, double, bool>()),                   \
        py::arg("particle_def"), py::arg("medium"), py::arg("energy_cuts"), py::arg("multiplier"), py::arg("lpm");

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

    // m_sub.doc() = R"pbdoc(
    //     pyPROPOSAL
    //     ----------
    //
    //     .. currentmodule:: component
    //
    //     .. autosummary::
    //         :toctree: _generate
    //
    //         Component
    //         Hydrogen
    //         Carbon
    //         Nitrogen
    //         Oxygen
    //         Sodium
    //         Magnesium
    //         Sulfur
    //         Argon
    //         Potassium
    //         Calcium
    //         Iron
    //         Copper
    //         Lead
    //         Uranium
    //         StandardRock
    //         FrejusRock
    // )pbdoc";

    py::class_<Components::Component, std::shared_ptr<Components::Component>>(m_sub, "Component")
        .def(py::init<std::string, double, double, double>(), py::arg("name"), py::arg("charge"), py::arg("atomic_num"), py::arg("atom_in_molecule"))
        .def("__str__", &py_print<Components::Component>)
        .def_property_readonly("name", &Components::Component::GetName)
        .def_property_readonly("nuclear_charge", &Components::Component::GetNucCharge)
        .def_property_readonly("atomic_number", &Components::Component::GetAtomicNum)
        .def_property_readonly("atoms_in_molecule", &Components::Component::GetAtomInMolecule)
        .def_property_readonly("log_constant", &Components::Component::GetLogConstant)
        .def_property_readonly("bprime", &Components::Component::GetBPrime)
        .def_property_readonly("average_nucleon_weight", &Components::Component::GetAverageNucleonWeight)
        .def_property_readonly("mn", &Components::Component::GetMN)
        .def_property_readonly("r0", &Components::Component::GetR0);

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
        .def_readonly("name", &ParticleDef::name)
        .def_readonly("mass", &ParticleDef::mass)
        .def_readonly("low", &ParticleDef::low)
        .def_readonly("charge", &ParticleDef::charge)
        .def_readonly("decay_table", &ParticleDef::decay_table)
        // .def_readonly("harc_component_table", &ParticleDef::hard_component_table)
        // .add_property("hard_component_table", make_function(&get_hard_component, return_internal_reference<>()))
        // //TODO(mario): shit Fri 2017/10/13
        ;

    py::class_<ParticleDef::Builder, std::shared_ptr<ParticleDef::Builder>>(m_sub, "ParticleDefBuilder")
        .def(py::init<>())
        .def("SetName", &ParticleDef::Builder::SetName)
        .def("SetMass", &ParticleDef::Builder::SetMass)
        .def("SetLow", &ParticleDef::Builder::SetLow)
        .def("SetLifetime", &ParticleDef::Builder::SetLifetime)
        .def("SetCharge", &ParticleDef::Builder::SetCharge)
        .def("SetDecayTable", &ParticleDef::Builder::SetDecayTable)
        .def("SetParticleDef", &ParticleDef::Builder::SetParticleDef)
        .def("build", &ParticleDef::Builder::build);


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
        .value("ContinuousEnergyLoss", DynamicData::ContinuousEnergyLoss);

    py::class_<DynamicData, std::shared_ptr<DynamicData> >(m_sub, "DynamicData")
        .def(py::init<DynamicData::Type>())
        .def(py::init<const DynamicData&>())
        .def("__str__", &py_print<DynamicData>)
        .def_property_readonly("id", &DynamicData::GetTypeId)
        .def_property("position", &DynamicData::GetPosition, &DynamicData::SetPosition)
        .def_property("direction", &DynamicData::GetDirection, &DynamicData::SetDirection)
        .def_property("energy", &DynamicData::GetEnergy, &DynamicData::SetEnergy)
        .def_property("parent_particle_energy", &DynamicData::GetParentParticleEnergy, &DynamicData::SetParentParticleEnergy)
        .def_property("time", &DynamicData::GetTime, &DynamicData::SetTime)
        .def_property("propagated_distance", &DynamicData::GetPropagatedDistance, &DynamicData::SetPropagatedDistance);

    py::class_<Particle, std::shared_ptr<Particle>, DynamicData>(m_sub, "Particle")
        .def(py::init<>())
        .def(py::init<const ParticleDef&>())
        .def(py::init<const Particle&>())
        .def("inject_state", &Particle::InjectState)
        .def_property_readonly("particle_def", &Particle::GetParticleDef)
        .def_property_readonly("decay_table", &Particle::GetDecayTable)
        .def_property("momentum", &Particle::GetMomentum, &Particle::SetMomentum)
        .def_property("entry_point", &Particle::GetEntryPoint, &Particle::SetEntryPoint)
        .def_property("entry_time", &Particle::GetEntryTime, &Particle::SetEntryTime)
        .def_property("entry_energy", &Particle::GetEntryEnergy, &Particle::SetEntryEnergy)
        .def_property("exit_point", &Particle::GetExitPoint, &Particle::SetExitPoint)
        .def_property("exit_time", &Particle::GetExitTime, &Particle::SetExitTime)
        .def_property("exit_energy", &Particle::GetExitEnergy, &Particle::SetExitEnergy)
        .def_property("closet_approach_point", &Particle::GetClosestApproachPoint, &Particle::SetClosestApproachPoint)
        .def_property("closet_approach_time", &Particle::GetClosestApproachTime, &Particle::SetClosestApproachTime)
        .def_property("closet_approach_energy", &Particle::GetClosestApproachEnergy, &Particle::SetClosestApproachEnergy)
        .def_property("e_lost", &Particle::GetElost, &Particle::SetElost);
}

void init_decay(py::module& m)
{
    py::module m_sub = m.def_submodule("decay");

    py::class_<DecayChannel, std::shared_ptr<DecayChannel>>(m_sub, "DecayChannel")
        .def("__str__", &py_print<DecayChannel>)
        .def("__eq__", &DecayChannel::operator==)
        .def("__ne__", &DecayChannel::operator!=)
        .def("decay", &DecayChannel::Decay, "Decay the given particle")
        .def_static("boost", (void (*)(Particle&, const Vector3D&, double)) &DecayChannel::Boost, "Boost the particle along a direction");

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

    py::enum_<GeometryFactory::Enum>(m_sub, "Shape")
        .value("Sphere", GeometryFactory::Sphere)
        .value("Box", GeometryFactory::Box)
        .value("Cylinder", GeometryFactory::Cylinder);

    py::class_<GeometryFactory::Definition, std::shared_ptr<GeometryFactory::Definition>>(m_sub, "GeometryDefinition")
        .def(py::init<>())
        .def_readwrite("shape", &GeometryFactory::Definition::shape)
        .def_readwrite("position", &GeometryFactory::Definition::position)
        .def_readwrite("inner_radius", &GeometryFactory::Definition::inner_radius)
        .def_readwrite("outer_radius", &GeometryFactory::Definition::radius)
        .def_readwrite("width", &GeometryFactory::Definition::width)
        .def_readwrite("height", &GeometryFactory::Definition::height)
        .def_readwrite("depth", &GeometryFactory::Definition::depth);

    py::class_<Geometry, std::shared_ptr<Geometry>>(m_sub, "Geometry")
        DEF_PY_PRINT(Geometry)
        .def("is_infront", &Geometry::IsInfront)
        .def("is_inside", &Geometry::IsInside)
        .def("is_behind", &Geometry::IsBehind)
        .def("distance_to_border", &Geometry::DistanceToBorder)
        .def("distance_to_closet_approach", &Geometry::DistanceToClosestApproach)
        .def_property_readonly("name", &Geometry::GetName)
        .def_property("position", &Geometry::GetPosition, &Geometry::SetPosition)
        .def_property("hirarchy", &Geometry::GetHirarchy, &Geometry::SetHirarchy);

    py::class_<Sphere, std::shared_ptr<Sphere>, Geometry>(m_sub, "Sphere")
        .def(py::init<>())
        .def(py::init<Vector3D, double, double>())
        .def(py::init<const Sphere&>())
        .def_property("inner_radius", &Sphere::GetInnerRadius, &Sphere::SetInnerRadius)
        .def_property("radius", &Sphere::GetRadius, &Sphere::SetRadius);

    py::class_<Box, std::shared_ptr<Box>, Geometry>(m_sub, "Box")
        .def(py::init<>())
        .def(py::init<Vector3D, double, double, double>())
        .def(py::init<const Box&>())
        .def_property("width", &Box::GetX, &Box::SetX)
        .def_property("height", &Box::GetY, &Box::SetY)
        .def_property("depth", &Box::GetZ, &Box::SetZ);

    py::class_<Cylinder, std::shared_ptr<Cylinder>, Geometry>(m_sub, "Cylinder")
        .def(py::init<>())
        .def(py::init<Vector3D, double, double, double>())
        .def(py::init<const Cylinder&>())
        .def_property("inner_radius", &Cylinder::GetInnerRadius, &Cylinder::SetInnerRadius)
        .def_property("radius", &Cylinder::GetRadius, &Cylinder::SetRadius)
        .def_property("height", &Cylinder::GetZ, &Cylinder::SetZ);
}

void init_parametrization(py::module& m)
{
    py::module m_sub = m.def_submodule("parametrization");

    // py::class_<Parametrization::IntegralLimits, std::shared_ptr<Parametrization::IntegralLimits>>("IntegralLimits")
    py::class_<Parametrization::IntegralLimits, std::shared_ptr<Parametrization::IntegralLimits>>(m_sub, "IntegralLimits")
        .def(py::init<>())
        .def_readwrite("v_max", &Parametrization::IntegralLimits::vMax)
        .def_readwrite("v_up", &Parametrization::IntegralLimits::vUp)
        .def_readwrite("v_min", &Parametrization::IntegralLimits::vMin);

    py::class_<Parametrization, std::shared_ptr<Parametrization>>(m_sub, "Parametrization")
        .def("__str__", &py_print<Parametrization>)
        .def("differential_crosssection", &Parametrization::DifferentialCrossSection)
        .def("dEdx_integrand", &Parametrization::FunctionToDEdxIntegral)
        .def("dE2dx_integrand", &Parametrization::FunctionToDE2dxIntegral)
        .def("dNdx_integrand", &Parametrization::FunctionToDNdxIntegral)
        .def("integral_limits", &Parametrization::GetIntegralLimits)
        .def_property_readonly("name", &Parametrization::GetName)
        .def_property_readonly("particle_def", &Parametrization::GetParticleDef)
        .def_property_readonly("medium", &Parametrization::GetMedium)
        .def_property_readonly("energy_cuts", &Parametrization::GetEnergyCuts)
        .def_property_readonly("multiplier", &Parametrization::GetMultiplier)
        .def_property_readonly("hash", &Parametrization::GetHash);

    // --------------------------------------------------------------------- //
    // Bremsstrahlung
    // --------------------------------------------------------------------- //

    py::module m_sub_brems = m_sub.def_submodule("bremsstrahlung");
    py::class_<Bremsstrahlung, std::shared_ptr<Bremsstrahlung>, Parametrization>(m_sub_brems, "Bremsstrahlung");

    BREMS_DEF(m_sub_brems, KelnerKokoulinPetrukhin)
    BREMS_DEF(m_sub_brems, PetrukhinShestakov)
    BREMS_DEF(m_sub_brems, CompleteScreening)
    BREMS_DEF(m_sub_brems, AndreevBezrukovBugaev)
    BREMS_DEF(m_sub_brems, SandrockSoedingreksoRhode)

    py::enum_<BremsstrahlungFactory::Enum>(m_sub_brems, "BremsParametrization")
        .value("PetrukhinShestakov", BremsstrahlungFactory::PetrukhinShestakov)
        .value("KelnerKokoulinPetrukhin", BremsstrahlungFactory::KelnerKokoulinPetrukhin)
        .value("CompleteScreening", BremsstrahlungFactory::CompleteScreening)
        .value("AndreevBezrukovBugaev", BremsstrahlungFactory::AndreevBezrukovBugaev);

    py::class_<BremsstrahlungFactory::Definition, std::shared_ptr<BremsstrahlungFactory::Definition> >(m_sub_brems, "BremsDefinition")
        .def(py::init<>())
        .def_readwrite("parametrization", &BremsstrahlungFactory::Definition::parametrization)
        .def_readwrite("lpm_effect", &BremsstrahlungFactory::Definition::lpm_effect)
        .def_readwrite("multiplier", &BremsstrahlungFactory::Definition::multiplier);

    // --------------------------------------------------------------------- //
    // Epair
    // --------------------------------------------------------------------- //

    py::module m_sub_epair = m_sub.def_submodule("pairproduction");
    py::class_<EpairProduction, std::shared_ptr<EpairProduction>, Parametrization>(m_sub_brems, "EpairProduction");

    py::class_<EpairProductionRhoIntegral, std::shared_ptr<EpairProductionRhoIntegral>, EpairProduction>(m_sub_epair, "EpairProductionRhoIntegral")
        .def("function_to_integral", &EpairProductionRhoIntegral::FunctionToIntegral);

    EPAIR_DEF(m_sub_epair, KelnerKokoulinPetrukhin)
    EPAIR_DEF(m_sub_epair, SandrockSoedingreksoRhode)

    EPAIR_INTERPOL_DEF(m_sub_epair, KelnerKokoulinPetrukhin)
    EPAIR_INTERPOL_DEF(m_sub_epair, SandrockSoedingreksoRhode)

    py::enum_<EpairProductionFactory::Enum>(m_sub_epair, "EpairParametrization")
        .value("PetrukhinShestakov", EpairProductionFactory::KelnerKokoulinPetrukhin)
        .value("KelnerKokoulinPetrukhin", EpairProductionFactory::SandrockSoedingreksoRhode);

    py::class_<EpairProductionFactory::Definition, std::shared_ptr<EpairProductionFactory::Definition> >(m_sub_epair, "EpairDefinition")
        .def(py::init<>())
        .def_readwrite("parametrization", &EpairProductionFactory::Definition::parametrization)
        .def_readwrite("lpm_effect", &EpairProductionFactory::Definition::lpm_effect)
        .def_readwrite("multiplier", &EpairProductionFactory::Definition::multiplier);

    // --------------------------------------------------------------------- //
    // Photo
    // --------------------------------------------------------------------- //

    py::module m_sub_photo = m_sub.def_submodule("photonuclear");
    py::class_<Photonuclear, std::shared_ptr<Photonuclear>, Parametrization>(m_sub_photo, "Photonuclear");

    // Shadow Effect
    py::class_<ShadowEffect, std::shared_ptr<ShadowEffect>>(m_sub_photo, "ShadowEffect")
        .def("calculate_shadow_effect", &ShadowEffect::CalculateShadowEffect)
        .def_property_readonly("name", &ShadowEffect::GetName);

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

    py::class_<PhotoRealPhotonAssumption, std::shared_ptr<PhotoRealPhotonAssumption>, Photonuclear>(m_sub_photo, "PhotoRealPhotonAssumption");
    py::class_<PhotoQ2Integral, std::shared_ptr<PhotoQ2Integral>, Photonuclear>(m_sub_photo, "PhotoQ2Integral");

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

    py::enum_<PhotonuclearFactory::Shadow>(m_sub_brems, "PhotoShadow")
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
             py::arg("multiplier"));

    py::class_<IonizationFactory::Definition, std::shared_ptr<IonizationFactory::Definition> >(m_sub_ioniz, "IonizationDefinition")
        .def(py::init<>())
        .def_readwrite("multiplier", &IonizationFactory::Definition::multiplier);
}

void init_crosssection(py::module& m)
{
    py::module m_sub = m.def_submodule("crosssection");

    py::class_<CrossSection, std::shared_ptr<CrossSection>>(m_sub, "CrossSection")
        DEF_PY_PRINT(CrossSection)
        .def("calculate_dEdx", &CrossSection::CalculatedEdx)
        .def("calculate_dE2dx", &CrossSection::CalculatedE2dx)
        .def("calculate_dNdx", (double (CrossSection::*)(double))&CrossSection::CalculatedNdx)
        .def("calculate_dNdx_rnd", (double (CrossSection::*)(double, double))&CrossSection::CalculatedNdx)
        .def("calculate_stochastic_loss", (double (CrossSection::*)(double, double, double))&CrossSection::CalculateStochasticLoss)
        .def_property_readonly("id", &CrossSection::GetTypeId)
        .def_property_readonly("parametrization", &CrossSection::GetParametrization);

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
        pyPROPOSAL
        ----------

        .. currentmodule:: pyPROPOSAL

        .. autosummary::
           :toctree: _generate

           Vector3D
           EnergyCutSettings
           InterpolationDef
           Utility
           UtilityDefinition
           ParticleLocation
           SectorDefinition
           Sector
           RandomGenerator
           Propagator
           PropagatorService
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

    py::class_<EnergyCutSettings, std::shared_ptr<EnergyCutSettings> >(m, "EnergyCutSettings")
        .def(py::init<>())
        .def(py::init<double, double>(), py::arg("ecut"), py::arg("vcut"))
        .def(py::init<const EnergyCutSettings&>())
        .def("__str__", &py_print<EnergyCutSettings>)
        .def_property("ecut", &EnergyCutSettings::GetEcut, &EnergyCutSettings::SetEcut)
        .def_property("vcut", &EnergyCutSettings::GetVcut, &EnergyCutSettings::SetVcut)
        .def("get_cut", &EnergyCutSettings::GetCut, "Return the lower from E*v = e");

    py::class_<InterpolationDef, std::shared_ptr<InterpolationDef>>(m, "InterpolationDef")
        .def(py::init<>())
        .def_readwrite("order_of_interpolation", &InterpolationDef::order_of_interpolation)
        .def_readwrite("path_to_tables", &InterpolationDef::path_to_tables)
        .def_readwrite("raw", &InterpolationDef::raw);

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
        .def_readwrite("ioniz_def", &Utility::Definition::ioniz_def);

    // --------------------------------------------------------------------- //
    // Sector Definition
    // --------------------------------------------------------------------- //

    // ----[ Location ]-------------------------------------- //

    py::enum_<Sector::ParticleLocation::Enum>(m, "ParticleLocation")
        .value("infront_detector", Sector::ParticleLocation::InfrontDetector)
        .value("inside_detector", Sector::ParticleLocation::InsideDetector)
        .value("behind_detector", Sector::ParticleLocation::BehindDetector);

    py::class_<Sector::Definition, std::shared_ptr<Sector::Definition> >(m, "SectorDefinition")
        .def(py::init<>())
        .def_readwrite("cut_settings", &Sector::Definition::cut_settings)
        .def_property("medium", &Sector::Definition::GetMedium, &Sector::Definition::SetMedium)
        .def_property("geometry", &Sector::Definition::GetGeometry, &Sector::Definition::SetGeometry)
        .def_readwrite("do_stochastic_loss_weighting", &Sector::Definition::do_stochastic_loss_weighting)
        .def_readwrite("stochastic_loss_weighting", &Sector::Definition::stochastic_loss_weighting)
        .def_readwrite("stopping_decay", &Sector::Definition::stopping_decay)
        .def_readwrite("do_continuous_randomization", &Sector::Definition::do_continuous_randomization)
        .def_readwrite("do_continuous_energy_loss_output", &Sector::Definition::do_continuous_energy_loss_output)
        .def_readwrite("do_exact_time_calculation", &Sector::Definition::do_exact_time_calculation)
        .def_readwrite("scattering_model", &Sector::Definition::scattering_model)
        .def_readwrite("particle_location", &Sector::Definition::location)
        .def_readwrite("crosssection_defs", &Sector::Definition::utility_def);

    // --------------------------------------------------------------------- //
    // Sector
    // --------------------------------------------------------------------- //

    py::class_<Sector, std::shared_ptr<Sector> >(m, "Sector")
        .def(py::init<Particle&, const Sector::Definition&>(), py::arg("particle"), py::arg("sector_definition"))
        .def(py::init<Particle&, const Sector::Definition&, const InterpolationDef&>(),
            py::arg("particle"), py::arg("sector_definition"), py::arg("interpolation_def"))
        .def("propagate", &Sector::Propagate, py::arg("distance"))
        .def("CalculateEnergyTillStochastic", &Sector::CalculateEnergyTillStochastic, py::arg("initial_energy"))
        .def_property_readonly("particle", &Sector::GetParticle, "Get the internal created particle to modify its properties");

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
        .def(py::init<const ParticleDef&, const std::vector<Sector::Definition>&, const Geometry&>(),
            py::arg("partcle_def"), py::arg("sector_defs"), py::arg("detector"))
        .def(py::init<const ParticleDef&, const std::vector<Sector::Definition>&, const Geometry&, const InterpolationDef&>(),
            py::arg("particle_def"), py::arg("sector_defs"), py::arg("detector"), py::arg("interpolation_def"))
        .def(py::init<const ParticleDef&, const std::vector<Sector::Definition>&, const Geometry&>(),
            py::arg("particle_def"), py::arg("sector_defs"), py::arg("detector"))
        .def(py::init<const ParticleDef&, const std::string&>(), py::arg("particle_def"), py::arg("config_file"))
        .def("propagate", &Propagator::Propagate, py::arg("max_distance_cm") = 1e20, py::return_value_policy::reference)
        .def_property_readonly("particle", &Propagator::GetParticle, "Get the internal created particle to modify its properties")
        .def_property_readonly("detector", &Propagator::GetDetector, "Get the detector geometry");

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
