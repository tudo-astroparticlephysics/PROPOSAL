#include "PROPOSAL/PROPOSAL.h"
#include "pyBindings.h"

#define BREMS_DEF(module, cls)                                           \
    py::class_<Brems##cls, std::shared_ptr<Brems##cls>, Bremsstrahlung>( \
        module, #cls)                                                    \
        .def(py::init<const ParticleDef&, const Medium&,                 \
                      const EnergyCutSettings&, double, bool>(),         \
             py::arg("particle_def"), py::arg("medium"),                 \
             py::arg("energy_cuts"), py::arg("multiplier"),              \
             py::arg("lpm_effect"));

#define PHOTO_REAL_DEF(module, cls, parent)                                    \
    py::class_<Photo##cls, std::shared_ptr<Photo##cls>, Photo##parent>(module, \
                                                                       #cls)   \
        .def(py::init<const ParticleDef&, const Medium&,                       \
                      const EnergyCutSettings&, double, bool>(),               \
             py::arg("particle_def"), py::arg("medium"),                       \
             py::arg("energy_cuts"), py::arg("multiplier"),                    \
             py::arg("add_pertubative"));

#define PHOTO_Q2_DEF(module, cls)                                              \
    py::class_<Photo##cls, std::shared_ptr<Photo##cls>, PhotoQ2Integral>(      \
        module, #cls)                                                          \
        .def(                                                                  \
            py::init<const ParticleDef&, const Medium&,                        \
                     const EnergyCutSettings&, double, const ShadowEffect&>(), \
            py::arg("particle_def"), py::arg("medium"),                        \
            py::arg("energy_cuts"), py::arg("multiplier"),                     \
            py::arg("shadow_effect"));

#define PHOTO_Q2_INTERPOL_DEF(module, cls)                                   \
    py::class_<PhotoQ2Interpolant<Photo##cls>,                               \
               std::shared_ptr<PhotoQ2Interpolant<Photo##cls>>, Photo##cls>( \
        module, #cls "Interpolant")                                          \
        .def(py::init<const ParticleDef&, const Medium&,                     \
                      const EnergyCutSettings&, double, const ShadowEffect&, \
                      InterpolationDef>(),                                   \
             py::arg("particle_def"), py::arg("medium"),                     \
             py::arg("energy_cuts"), py::arg("multiplier"),                  \
             py::arg("shadow_effect"), py::arg("interpolation_def"));

#define EPAIR_DEF(module, cls)                                   \
    py::class_<Epair##cls, std::shared_ptr<Epair##cls>,          \
               EpairProductionRhoIntegral>(module, #cls)         \
        .def(py::init<const ParticleDef&, const Medium&,         \
                      const EnergyCutSettings&, double, bool>(), \
             py::arg("particle_def"), py::arg("medium"),         \
             py::arg("energy_cuts"), py::arg("multiplier"),      \
             py::arg("lpm_effect"));

#define EPAIR_INTERPOL_DEF(module, cls)                                    \
    py::class_<EpairProductionRhoInterpolant<Epair##cls>,                  \
               std::shared_ptr<EpairProductionRhoInterpolant<Epair##cls>>, \
               Epair##cls>(module, #cls "Interpolant")                     \
        .def(py::init<const ParticleDef&, const Medium&,                   \
                      const EnergyCutSettings&, double, bool,              \
                      InterpolationDef>(),                                 \
             py::arg("particle_def"), py::arg("medium"),                   \
             py::arg("energy_cuts"), py::arg("multiplier"),                \
             py::arg("lpm_effect"), py::arg("interpolation_def"));

#define MUPAIR_DEF(module, cls)                            \
    py::class_<Mupair##cls, std::shared_ptr<Mupair##cls>,  \
               MupairProductionRhoIntegral>(module, #cls)  \
        .def(py::init<const ParticleDef&, const Medium&,   \
                      const EnergyCutSettings&, double>(), \
             py::arg("particle_def"), py::arg("medium"),   \
             py::arg("energy_cuts"), py::arg("multiplier"));

#define MUPAIR_INTERPOL_DEF(module, cls)                                     \
    py::class_<MupairProductionRhoInterpolant<Mupair##cls>,                  \
               std::shared_ptr<MupairProductionRhoInterpolant<Mupair##cls>>, \
               Mupair##cls>(module, #cls "Interpolant")                      \
        .def(py::init<const ParticleDef&, const Medium&,                     \
                      const EnergyCutSettings&, double, InterpolationDef>(), \
             py::arg("particle_def"), py::arg("medium"),                     \
             py::arg("energy_cuts"), py::arg("multiplier"),                  \
             py::arg("interpolation_def"));

namespace py = pybind11;
using namespace PROPOSAL;

void init_parametrization(py::module& m) {
    py::module m_sub = m.def_submodule("parametrization");

    // py::class_<Parametrization::IntegralLimits,
    // std::shared_ptr<Parametrization::IntegralLimits>>("IntegralLimits")
    py::class_<Parametrization::IntegralLimits,
               std::shared_ptr<Parametrization::IntegralLimits>>(
        m_sub, "IntegralLimits")
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

    py::class_<Parametrization, std::shared_ptr<Parametrization>>(
        m_sub, "Parametrization",
        R"pbdoc(
            Parametrization objects provide the theoretical input for physical cross section used in PROPOSAL, whereas :meth:`~pyPROPOSAL.crosssection.CrossSection`
            provides the numerical methods to process the parametrization. 

            For each physical process in PROPOSAL there are several different parametrizations available, so the user can check how the theoretical input influences
            the simulation.
            )pbdoc")
        .def("__str__", &py_print<Parametrization>)
        .def("differential_crosssection",
             &Parametrization::DifferentialCrossSection, py::arg("energy"),
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

                )pbdoc")
        .def("dEdx_integrand", &Parametrization::FunctionToDEdxIntegral,
             py::arg("energy"), py::arg("v"),
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

                )pbdoc")
        .def("dE2dx_integrand", &Parametrization::FunctionToDE2dxIntegral,
             py::arg("energy"), py::arg("v"),
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

                )pbdoc")
        .def("dNdx_integrand", &Parametrization::FunctionToDNdxIntegral,
             py::arg("energy"), py::arg("v"),
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

                )pbdoc")
        .def("integral_limits", &Parametrization::GetIntegralLimits,
             py::arg("energy"),
             R"pbdoc(
            Returns:
                integral_limits (:meth:`~pyPROPOSAL.parametrization.IntegralLimits`): returns the integral limits for the given energy and the
                current parametrization
                )pbdoc")
        .def_property_readonly("name", &Parametrization::GetName,
                               R"pbdoc(

            Get name of current parametrization

                )pbdoc")
        .def_property_readonly("particle_def", &Parametrization::GetParticleDef,
                               R"pbdoc(

            Get :meth:`~pyPROPOSAL.particle.ParticleDef` used by the parametrization

                )pbdoc")
        .def_property_readonly("medium", &Parametrization::GetMedium,
                               R"pbdoc( 

            Get :meth:`~pyPROPOSAL.medium` used by the parametrization

                )pbdoc")
        .def_property_readonly("energy_cuts", &Parametrization::GetEnergyCuts,
                               R"pbdoc( 

            Get :meth:`~pyPROPOSAL.EnergyCutSettings` defined by the user for the parametrization

                )pbdoc")
        .def_property_readonly("multiplier", &Parametrization::GetMultiplier,
                               R"pbdoc( 

            Get multiplier used for the parametrization

                )pbdoc")
        .def_property_readonly("hash", &Parametrization::GetHash,
                               R"pbdoc( 

            Get internal hash corresponding to the current parametrization

                )pbdoc");

    // ---------------------------------------------------------------------
    // // Bremsstrahlung
    // ---------------------------------------------------------------------
    // //

    py::module m_sub_brems = m_sub.def_submodule("bremsstrahlung");
    py::class_<Bremsstrahlung, std::shared_ptr<Bremsstrahlung>,
               Parametrization>(m_sub_brems, "Bremsstrahlung",
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
                )pbdoc");

    BREMS_DEF(m_sub_brems, KelnerKokoulinPetrukhin)
    BREMS_DEF(m_sub_brems, PetrukhinShestakov)
    BREMS_DEF(m_sub_brems, CompleteScreening)
    BREMS_DEF(m_sub_brems, AndreevBezrukovBugaev)
    BREMS_DEF(m_sub_brems, SandrockSoedingreksoRhode)

    py::enum_<BremsstrahlungFactory::Enum>(m_sub_brems, "BremsParametrization")
        .value("PetrukhinShestakov", BremsstrahlungFactory::PetrukhinShestakov)
        .value("KelnerKokoulinPetrukhin",
               BremsstrahlungFactory::KelnerKokoulinPetrukhin)
        .value("CompleteScreening", BremsstrahlungFactory::CompleteScreening)
        .value("AndreevBezrukovBugaev",
               BremsstrahlungFactory::AndreevBezrukovBugaev)
        .value("SandrockSoedingreksoRhode",
               BremsstrahlungFactory::SandrockSoedingreksoRhode);

    py::class_<BremsstrahlungFactory,
               std::unique_ptr<BremsstrahlungFactory, py::nodelete>>(
        m_sub_brems, "BremsFactory")
        .def("get_enum_from_str", &BremsstrahlungFactory::GetEnumFromString,
             py::arg("parametrization_str"))
        .def("create_bremsstrahlung",
             (CrossSection *
              (BremsstrahlungFactory::*)(const ParticleDef&, const Medium&,
                                         const EnergyCutSettings&,
                                         const BremsstrahlungFactory::
                                             Definition&)const) &
                 BremsstrahlungFactory::CreateBremsstrahlung,
             py::arg("particle_def"), py::arg("medium"), py::arg("ecuts"),
             py::arg("brems_def"))
        .def("create_bremsstrahlung_interpol",
             (CrossSection *
              (BremsstrahlungFactory::*)(const ParticleDef&, const Medium&,
                                         const EnergyCutSettings&,
                                         const BremsstrahlungFactory::
                                             Definition&,
                                         InterpolationDef) const) &
                 BremsstrahlungFactory::CreateBremsstrahlung,
             py::arg("particle_def"), py::arg("medium"), py::arg("ecuts"),
             py::arg("brems_def"), py::arg("interpolation_def"))
        .def_static("get", &BremsstrahlungFactory::Get,
                    py::return_value_policy::reference);

    py::class_<BremsstrahlungFactory::Definition,
               std::shared_ptr<BremsstrahlungFactory::Definition>>(
        m_sub_brems, "BremsDefinition")
        .def(py::init<>())
        .def_readwrite("parametrization",
                       &BremsstrahlungFactory::Definition::parametrization)
        .def_readwrite("lpm_effect",
                       &BremsstrahlungFactory::Definition::lpm_effect)
        .def_readwrite("multiplier",
                       &BremsstrahlungFactory::Definition::multiplier);

    // ---------------------------------------------------------------------
    // // Epair
    // ---------------------------------------------------------------------
    // //

    py::module m_sub_epair = m_sub.def_submodule("pairproduction");
    py::class_<EpairProduction, std::shared_ptr<EpairProduction>,
               Parametrization>(m_sub_epair, "EpairProduction",
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

    py::class_<EpairProductionRhoIntegral,
               std::shared_ptr<EpairProductionRhoIntegral>, EpairProduction>(
        m_sub_epair, "EpairProductionRhoIntegral")
        .def("function_to_integral",
             &EpairProductionRhoIntegral::FunctionToIntegral);

    EPAIR_DEF(m_sub_epair, KelnerKokoulinPetrukhin)
    EPAIR_DEF(m_sub_epair, SandrockSoedingreksoRhode)

    EPAIR_INTERPOL_DEF(m_sub_epair, KelnerKokoulinPetrukhin)
    EPAIR_INTERPOL_DEF(m_sub_epair, SandrockSoedingreksoRhode)

    py::enum_<EpairProductionFactory::Enum>(m_sub_epair, "EpairParametrization")
        .value("KelnerKokoulinPetrukhin",
               EpairProductionFactory::KelnerKokoulinPetrukhin)
        .value("SandrockSoedingreksoRhode",
               EpairProductionFactory::SandrockSoedingreksoRhode);

    py::class_<EpairProductionFactory,
               std::unique_ptr<EpairProductionFactory, py::nodelete>>(
        m_sub_epair, "EpairFactory")
        .def("get_enum_from_str", &EpairProductionFactory::GetEnumFromString,
             py::arg("parametrization_str"))
        .def("create_pairproduction",
             (CrossSection *
              (EpairProductionFactory::*)(const ParticleDef&, const Medium&,
                                          const EnergyCutSettings&,
                                          const EpairProductionFactory::
                                              Definition&)const) &
                 EpairProductionFactory::CreateEpairProduction,
             py::arg("particle_def"), py::arg("medium"), py::arg("ecuts"),
             py::arg("brems_def"))
        .def("create_pairproduction_interpol",
             (CrossSection *
              (EpairProductionFactory::*)(const ParticleDef&, const Medium&,
                                          const EnergyCutSettings&,
                                          const EpairProductionFactory::
                                              Definition&,
                                          InterpolationDef) const) &
                 EpairProductionFactory::CreateEpairProduction,
             py::arg("particle_def"), py::arg("medium"), py::arg("ecuts"),
             py::arg("brems_def"), py::arg("interpolation_def"))
        .def_static("get", &EpairProductionFactory::Get,
                    py::return_value_policy::reference);

    py::class_<EpairProductionFactory::Definition,
               std::shared_ptr<EpairProductionFactory::Definition>>(
        m_sub_epair, "EpairDefinition")
        .def(py::init<>())
        .def_readwrite("parametrization",
                       &EpairProductionFactory::Definition::parametrization)
        .def_readwrite("lpm_effect",
                       &EpairProductionFactory::Definition::lpm_effect)
        .def_readwrite("multiplier",
                       &EpairProductionFactory::Definition::multiplier);

    // ---------------------------------------------------------------------
    // // Mupair
    // ---------------------------------------------------------------------
    // //

    py::module m_sub_mupair = m_sub.def_submodule("mupairproduction");
    py::class_<MupairProduction, std::shared_ptr<MupairProduction>,
               Parametrization>(m_sub_mupair, "MupairProduction",
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
        .def("Calculaterho", &MupairProduction::Calculaterho, py::arg("energy"),
             py::arg("v"), py::arg("rnd1"), py::arg("rnd2"),
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

                )pbdoc");

    py::class_<MupairProductionRhoIntegral,
               std::shared_ptr<MupairProductionRhoIntegral>, MupairProduction>(
        m_sub_mupair, "MupairProductionRhoIntegral")
        .def("function_to_integral",
             &MupairProductionRhoIntegral::FunctionToIntegral);

    MUPAIR_DEF(m_sub_mupair, KelnerKokoulinPetrukhin)

    MUPAIR_INTERPOL_DEF(m_sub_mupair, KelnerKokoulinPetrukhin)

    py::enum_<MupairProductionFactory::Enum>(m_sub_mupair,
                                             "MupairParametrization")
        .value("KelnerKokoulinPetrukhin",
               MupairProductionFactory::KelnerKokoulinPetrukhin);

    py::class_<MupairProductionFactory::Definition,
               std::shared_ptr<MupairProductionFactory::Definition>>(
        m_sub_mupair, "MupairDefinition")
        .def(py::init<>())
        .def_readwrite("parametrization",
                       &MupairProductionFactory::Definition::parametrization)
        .def_readwrite("mupair_enable",
                       &MupairProductionFactory::Definition::mupair_enable)
        .def_readwrite("multiplier",
                       &MupairProductionFactory::Definition::multiplier)
        .def_readwrite("particle_output",
                       &MupairProductionFactory::Definition::particle_output);

    // ---------------------------------------------------------------------
    // // Weak Interaction
    // ---------------------------------------------------------------------
    // //

    py::module m_sub_weak = m_sub.def_submodule("weakinteraction");
    py::class_<WeakInteraction, std::shared_ptr<WeakInteraction>,
               Parametrization>(m_sub_weak, "WeakInteraction",
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

    py::class_<WeakCooperSarkarMertsch,
               std::shared_ptr<WeakCooperSarkarMertsch>, WeakInteraction>(
        m_sub_weak, "CooperSarkarMertsch")
        .def(py::init<const ParticleDef&, const Medium&, double>(),
             py::arg("particle_def"), py::arg("medium"), py::arg("multiplier"));

    py::enum_<WeakInteractionFactory::Enum>(m_sub_weak, "WeakParametrization")
        .value("CooperSarkarMertsch",
               WeakInteractionFactory::CooperSarkarMertsch);

    py::class_<WeakInteractionFactory::Definition,
               std::shared_ptr<WeakInteractionFactory::Definition>>(
        m_sub_weak, "WeakDefinition")
        .def(py::init<>())
        .def_readwrite("parametrization",
                       &WeakInteractionFactory::Definition::parametrization)
        .def_readwrite("weak_enable",
                       &WeakInteractionFactory::Definition::weak_enable)
        .def_readwrite("multiplier",
                       &WeakInteractionFactory::Definition::multiplier);

    // ---------------------------------------------------------------------
    // // Photo
    // ---------------------------------------------------------------------
    // //

    py::module m_sub_photo = m_sub.def_submodule("photonuclear");
    py::class_<Photonuclear, std::shared_ptr<Photonuclear>, Parametrization>(
        m_sub_photo, "Photonuclear");

    // Shadow Effect
    py::class_<ShadowEffect, std::shared_ptr<ShadowEffect>>(m_sub_photo,
                                                            "ShadowEffect",
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

    py::class_<ShadowDuttaRenoSarcevicSeckel,
               std::shared_ptr<ShadowDuttaRenoSarcevicSeckel>, ShadowEffect>(
        m_sub_photo, "ShadowDuttaRenoSarcevicSeckel")
        .def(py::init<>());

    py::class_<ShadowButkevichMikhailov,
               std::shared_ptr<ShadowButkevichMikhailov>, ShadowEffect>(
        m_sub_photo, "ShadowButkevichMikhailov")
        .def(py::init<>());

    // Real Photon
    py::class_<RealPhoton, std::shared_ptr<RealPhoton>>(m_sub_photo,
                                                        "RealPhoton")
        .def("calculate_hard_component", &RealPhoton::CalculateHardComponent)
        .def_property_readonly("name", &RealPhoton::GetName);

    py::class_<SoftComponent, std::shared_ptr<SoftComponent>, RealPhoton>(
        m_sub_photo, "SoftComponent")
        .def(py::init<>());
    py::class_<HardComponent, std::shared_ptr<HardComponent>>(m_sub_photo,
                                                              "HardComponent")
        .def(py::init<const ParticleDef&>(), py::arg("particle_def"));

    py::class_<PhotoRealPhotonAssumption,
               std::shared_ptr<PhotoRealPhotonAssumption>, Photonuclear>(
        m_sub_photo, "PhotoRealPhotonAssumption",
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
    py::class_<PhotoQ2Integral, std::shared_ptr<PhotoQ2Integral>, Photonuclear>(
        m_sub_photo, "PhotoQ2Integral",
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
    PHOTO_REAL_DEF(m_sub_photo, Kokoulin,
                   BezrukovBugaev)  // Kokoulin derives from BezrukovBugaev

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
        .value("AbramowiczLevinLevyMaor91",
               PhotonuclearFactory::AbramowiczLevinLevyMaor91)
        .value("AbramowiczLevinLevyMaor97",
               PhotonuclearFactory::AbramowiczLevinLevyMaor97)
        .value("ButkevichMikhailov", PhotonuclearFactory::ButkevichMikhailov)
        .value("RenoSarcevicSu", PhotonuclearFactory::RenoSarcevicSu);

    py::enum_<PhotonuclearFactory::Shadow>(m_sub_photo, "PhotoShadow")
        .value("DuttaRenoSarcevicSeckel",
               PhotonuclearFactory::ShadowDuttaRenoSarcevicSeckel)
        .value("ButkevichMikhailov",
               PhotonuclearFactory::ShadowButkevichMikhailov);

    py::class_<PhotonuclearFactory,
               std::unique_ptr<PhotonuclearFactory, py::nodelete>>(
        m_sub_photo, "PhotoFactory")
        .def("get_enum_from_str", &PhotonuclearFactory::GetEnumFromString,
             py::arg("parametrization_str"))
        .def("create_photonuclear",
             (CrossSection *
              (PhotonuclearFactory::*)(const ParticleDef&, const Medium&,
                                       const EnergyCutSettings&,
                                       const PhotonuclearFactory::Definition&)
                  const) &
                 PhotonuclearFactory::CreatePhotonuclear,
             py::arg("particle_def"), py::arg("medium"), py::arg("ecuts"),
             py::arg("photo_def"))
        .def("create_photonuclear_interpol",
             (CrossSection *
              (PhotonuclearFactory::*)(const ParticleDef&, const Medium&,
                                       const EnergyCutSettings&,
                                       const PhotonuclearFactory::Definition&,
                                       InterpolationDef) const) &
                 PhotonuclearFactory::CreatePhotonuclear,
             py::arg("particle_def"), py::arg("medium"), py::arg("ecuts"),
             py::arg("photo_def"), py::arg("interpolation_def"))
        .def_static("get", &PhotonuclearFactory::Get,
                    py::return_value_policy::reference);

    py::class_<PhotonuclearFactory::Definition,
               std::shared_ptr<PhotonuclearFactory::Definition>>(
        m_sub_photo, "PhotoDefinition")
        .def(py::init<>())
        .def_readwrite("parametrization",
                       &PhotonuclearFactory::Definition::parametrization)
        .def_readwrite("shadow", &PhotonuclearFactory::Definition::shadow)
        .def_readwrite("hard_component",
                       &PhotonuclearFactory::Definition::hard_component)
        .def_readwrite("multiplier",
                       &PhotonuclearFactory::Definition::multiplier);

    // ---------------------------------------------------------------------
    // // Ionization
    // ---------------------------------------------------------------------
    // //

    py::module m_sub_ioniz = m_sub.def_submodule("ionization");

    py::class_<Ionization, std::shared_ptr<Ionization>, Parametrization>(
        m_sub_ioniz, "Ionization")
        .def(py::init<const ParticleDef&, const Medium&,
                      const EnergyCutSettings&, double>(),
             py::arg("particle_def"), py::arg("medium"), py::arg("energy_cuts"),
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

    py::class_<IonizationFactory::Definition,
               std::shared_ptr<IonizationFactory::Definition>>(
        m_sub_ioniz, "IonizationDefinition")
        .def(py::init<>())
        .def_readwrite("multiplier",
                       &IonizationFactory::Definition::multiplier);
}
