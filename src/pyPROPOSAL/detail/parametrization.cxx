#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/crosssection/parametrization/Annihilation.h"
#include "PROPOSAL/crosssection/parametrization/Bremsstrahlung.h"
#include "PROPOSAL/crosssection/parametrization/Compton.h"
#include "PROPOSAL/crosssection/parametrization/EpairProduction.h"
#include "PROPOSAL/crosssection/parametrization/Ionization.h"
#include "PROPOSAL/crosssection/parametrization/MupairProduction.h"
#include "PROPOSAL/crosssection/parametrization/PhotoPairProduction.h"
#include "PROPOSAL/crosssection/parametrization/PhotoQ2Integration.h"
#include "PROPOSAL/crosssection/parametrization/PhotoRealPhotonAssumption.h"
#include "PROPOSAL/crosssection/parametrization/Photonuclear.h"
#include "PROPOSAL/crosssection/parametrization/WeakInteraction.h"
#include "pyPROPOSAL/pyBindings.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/medium/Medium.h"

#define BREMS_DEF(module, cls)                                                 \
    py::class_<crosssection::Brems##cls,                                       \
        std::shared_ptr<crosssection::Brems##cls>,                             \
        crosssection::Bremsstrahlung>(module, #cls)                            \
        .def(py::init<bool>(), py::arg("lpm") = false)                         \
        .def(py::init<bool, const ParticleDef&, const Medium&, double>(),      \
            py::arg("lpm"), py::arg("particle_def"), py::arg("medium"),        \
            py::arg("density_correction") = 1.0);

#define PHOTO_REAL_DEF(module, cls, parent)                                    \
    py::class_<crosssection::Photo##cls,                                       \
        std::shared_ptr<crosssection::Photo##cls>,                             \
        crosssection::Photo##parent>(module, #cls)                             \
        .def(py::init<bool>(), py::arg("add_pertubative"));

#define PHOTO_Q2_DEF(module, cls)                                              \
    py::class_<crosssection::Photo##cls,                                       \
        std::shared_ptr<crosssection::Photo##cls>,                             \
        crosssection::PhotoQ2Integral>(module, #cls)                           \
        .def(py::init<std::shared_ptr<crosssection::ShadowEffect>>(),          \
            py::arg("shadow_effect"));

#define EPAIR_DEF(module, cls)                                                 \
    py::class_<crosssection::Epair##cls,                                       \
        std::shared_ptr<crosssection::Epair##cls>,                             \
        crosssection::EpairProductionRhoIntegral>(module, #cls)                \
        .def(py::init<bool>(), py::arg("lpm") = false)                  \
        .def(py::init<bool, const ParticleDef&, const Medium&, double>(),      \
            py::arg("lpm"), py::arg("particle_def"), py::arg("medium"),        \
            py::arg("density_correction") = 1.0);

#define MUPAIR_DEF(module, cls)                                                \
    py::class_<crosssection::Mupair##cls,                                      \
        std::shared_ptr<crosssection::Mupair##cls>,                            \
        crosssection::MupairProductionRhoIntegral>(module, #cls)               \
        .def(py::init<>());

namespace py = pybind11;
using namespace PROPOSAL;

void init_parametrization(py::module& m)
{
    py::module m_sub = m.def_submodule("parametrization");

    auto param_docstring_class = R"pbdoc(
            Parametrization objects provide the theoretical input for physical
            cross section used in PROPOSAL, whereas
            :meth:`~proposal.crosssection.CrossSection` provides the numerical
            methods to process the parametrization.

            For each physical process in PROPOSAL there are several different
            parametrizations available, so the user can check how the
            theoretical input influences the simulation.
    )pbdoc";

    auto param_docstring_diff_cross = R"pbdoc(
            Calculate the value

            .. math::

                \frac{d\sigma}{dv}(E),

            e.g. the differential crosssection in v for the current
            parametrization.  If the parametrization is given in a
            double-differential-crosssection the crosssection will be
            integrated over all values but v.

            Args:
                energy (float): energy in MeV
                v (float): fraction of energy that is transfered from the
                    primary particle to secondary particles (or other forms of
                    energy losses)

            Return:
                differential_crosssection (float): returns the differential
                    crosssection in v
    )pbdoc";

    auto param_docstring_dEdx_integrand = R"pbdoc(
            Calculate the value

            .. math::

                v \cdot \frac{d\sigma}{dv}(E),

            e.g. the differential crosssection in v for the current
            parametrization multipied by v.  If the parametrization is given in
            a double-differential-crosssection the crosssection will be
            integrated over all values but v.

            Args:
                energy (float): energy in MeV
                v (float): fraction of energy that is transfered from the
                    primary particle to secondary particles (or other forms of
                    energy losses)

            Return:
                dEdx_integrand (float): returns the differential crosssection in
                    v multiplied by v

            PROPOSAL uses this function internally, for example to calcuate
                :math:`\langle\frac{dE}{dx}\rangle`
    )pbdoc";

    auto param_docstring_dE2dx_integrand = R"pbdoc(
            Calculate the value

            .. math::

                v^2 \cdot \frac{d\sigma}{dv}(E),

            e.g. the differential crosssection in v for the current
            parametrization multipied by :math:`v^2`.  If the parametrization is
            given in a double-differential-crosssection the crosssection will be
            integrated over all values but v.

            Args:
                energy (float): energy in MeV
                v (float): fraction of energy that is transfered from the
                    primary particle to secondary particles (or other forms of
                    energy losses)

            Return:
                dE2dx_integrand (float): returns the differential crosssection in
                    v multiplied by :math:`v^2`

            PROPOSAL uses this function internally, for example in the
            calculation of the countinous randomization.
    )pbdoc";

    auto param_docstring_kinematic_limits = R"pbdoc(
            Get internal hash corresponding to the current parametrization
    )pbdoc";

    py::class_<crosssection::Parametrization<Medium>,
        std::shared_ptr<crosssection::Parametrization<Medium>>>(
        m_sub, "ParametrizationForMedium", param_docstring_class)
        .def("differential_crosssection",
            py::overload_cast<ParticleDef const&, Medium const&,
                double, double>(
                &crosssection::Parametrization<Medium>::DifferentialCrossSection,
                py::const_),
            py::arg("particle_def"), py::arg("medium"), py::arg("energy"),
            py::arg("v"), param_docstring_diff_cross)
        .def("dEdx_integrand",
            &crosssection::Parametrization<Medium>::FunctionToDEdxIntegral,
            py::arg("particle_def"), py::arg("medium"), py::arg("energy"),
            py::arg("v"), param_docstring_dEdx_integrand)
        .def("dE2dx_integrand",
            &crosssection::Parametrization<Medium>::FunctionToDE2dxIntegral,
            py::arg("particle_def"), py::arg("medium"), py::arg("energy"),
            py::arg("v"), param_docstring_dE2dx_integrand)
        .def("kinematic_limits",
            &crosssection::Parametrization<Medium>::GetKinematicLimits,
            py::arg("particle_def"), py::arg("medium"), py::arg("energy"))
        .def_property_readonly("hash", &crosssection::Parametrization<Medium>::GetHash,
            param_docstring_kinematic_limits);

    py::class_<crosssection::Parametrization<Component>,
            std::shared_ptr<crosssection::Parametrization<Component>>>(
            m_sub, "ParametrizationForComponent", param_docstring_class)
            .def("differential_crosssection",
                 py::overload_cast<ParticleDef const&, Component const&,
                         double, double>(
                         &crosssection::Parametrization<Component>::DifferentialCrossSection,
                         py::const_),
                 py::arg("particle_def"), py::arg("component"), py::arg("energy"),
                 py::arg("v"), param_docstring_diff_cross)
            .def("dEdx_integrand",
                 &crosssection::Parametrization<Component>::FunctionToDEdxIntegral,
                 py::arg("particle_def"), py::arg("component"), py::arg("energy"),
                 py::arg("v"), param_docstring_dEdx_integrand)
            .def("dE2dx_integrand",
                 &crosssection::Parametrization<Component>::FunctionToDE2dxIntegral,
                 py::arg("particle_def"), py::arg("component"), py::arg("energy"),
                 py::arg("v"), param_docstring_dE2dx_integrand)
            .def("kinematic_limits",
                 &crosssection::Parametrization<Component>::GetKinematicLimits,
                 py::arg("particle_def"), py::arg("component"), py::arg("energy"))
            .def_property_readonly("hash", &crosssection::Parametrization<Component>::GetHash,
                                   param_docstring_kinematic_limits);

    // ---------------------------------------------------------------------
    // // Bremsstrahlung
    // ---------------------------------------------------------------------
    // //

    py::module m_sub_brems = m_sub.def_submodule("bremsstrahlung");
    py::class_<crosssection::Bremsstrahlung,
        std::shared_ptr<crosssection::Bremsstrahlung>,
        crosssection::Parametrization<Component>>(m_sub_brems, "Bremsstrahlung",
        R"pbdoc(

            Virtual class for the Bremsstrahlung parametrizations. They can be initialized by using one of the given parametrizations with the following parameters

            Args:
                particle_def (:meth:`~proposal.particle.ParticleDef`): includes all static particle information for the parametrization such as mass, charge, etc.
                medium (:meth:`~proposal.medium`): includes all medium information for the parametrization such as densities or nucleon charges
                energy_cuts (:meth:`~proposal.EnergyCutSettings`): energy cut setting for the parametrization
                multiplier (double): Use a multiplicative factor for the differential crosssection. Can be used for testing or other studies
                lpm (bool): Enable or disable the corrections due to the Ter-Mikaelian and Landau-Pomeranchuk effect.

            The following parametrizations are currently implemented:

            * KelnerKokoulinPetrukhin

            * PetrukhinShestakov

            * CompleteScreening

            * AndreevBezrukovBugaev

            * SandrockSoedingreksoRhode

            Example:
                To create a bremsstrahlung parametrization

                >>> mu = proposal.particle.MuMinusDef.get()
                >>> medium = proposal.medium.StandardRock(1.0)
                >>> cuts = proposal.EnergyCutSettings(-1, -1)
                >>> param = proposal.parametrization.bremsstrahlung.SandrockSoedingreksoRhode(mu, medium, cuts, 1.0, False)
                )pbdoc");

    BREMS_DEF(m_sub_brems, KelnerKokoulinPetrukhin)
    BREMS_DEF(m_sub_brems, PetrukhinShestakov)
    BREMS_DEF(m_sub_brems, CompleteScreening)
    BREMS_DEF(m_sub_brems, AndreevBezrukovBugaev)
    BREMS_DEF(m_sub_brems, SandrockSoedingreksoRhode)
    BREMS_DEF(m_sub_brems, ElectronScreening)

    py::class_<crosssection::BremsLPM, std::shared_ptr<crosssection::BremsLPM>>(
        m_sub_brems, "brems_lpm")
        .def(py::init<const ParticleDef&, const Medium&,
                 const crosssection::Bremsstrahlung&, double>(),
            py::arg("particle_def"), py::arg("medium"),
            py::arg("bremsstrahlung"), py::arg("density_correction") = 1.0)
        .def("supression_factor", &crosssection::BremsLPM::suppression_factor,
            py::arg("energy"), py::arg("v"), py::arg("component"));

    // ---------------------------------------------------------------------
    // // Epair
    // ---------------------------------------------------------------------
    // //

    py::module m_sub_epair = m_sub.def_submodule("pairproduction");
    py::class_<crosssection::EpairProduction,
        std::shared_ptr<crosssection::EpairProduction>,
        crosssection::Parametrization<Component>>(m_sub_epair, "EpairProduction",
        R"pbdoc(

            Virtual class for the electron pair production parametrizations. They can be initialized by using one of the given parametrizations with the following parameters

            Args:
                particle_def (:meth:`~proposal.particle.ParticleDef`): includes all static particle information for the parametrization such as mass, charge, etc.
                medium (:meth:`~proposal.medium`): includes all medium information for the parametrization such as densities or nucleon charges
                energy_cuts (:meth:`~proposal.EnergyCutSettings`): energy cut setting for the parametrization
                multiplier (double): Use a multiplicative factor for the differential crosssection. Can be used for testing or other studies
                lpm (bool): Enable or disable the corrections due to the Ter-Mikaelian and Landau-Pomeranchuk effect.
                interpolation_def (:meth:`~proposal.InterpolationDef`): Only needed by Interpolant parametrizations. Includes settings for the interpolation

            Since the differential cross section is given in :math:`\rho` as well, an intergration over this parameter is needed.
            When using the interpolation_def parameter, this integration is saved in interpolation tables (improving the performance of the calculation with neglible decline in accuracy).

            The following parametrizations are currently implemented:

            * KelnerKokoulinPetrukhin

            * SandrockSoedingreksoRhode

            * ForElectronPositron

            * KelnerKokoulinPetrukhinInterpolant

            * SandrockSoedingreksoRhodeInterpolant

            * ForElectronPositronInterpolant

            Example:
                To create a electron pair production parametrization

                >>> mu = proposal.particle.MuMinusDef.get()
                >>> medium = proposal.medium.StandardRock(1.0)
                >>> cuts = proposal.EnergyCutSettings(-1, -1)
                >>> param = proposal.parametrization.pairproduction.SandrockSoedingreksoRhode(mu, medium, cuts, 1.0, False)
                )pbdoc");

    py::class_<crosssection::EpairProductionRhoIntegral,
        std::shared_ptr<crosssection::EpairProductionRhoIntegral>,
        crosssection::EpairProduction>(
        m_sub_epair, "EpairProductionRhoIntegral")
        .def("function_to_integral",
            &crosssection::EpairProductionRhoIntegral::FunctionToIntegral,
             py::arg("particle_def"), py::arg("component"), py::arg("energy"),
             py::arg("v"), py::arg("rho"));

    EPAIR_DEF(m_sub_epair, KelnerKokoulinPetrukhin)
    EPAIR_DEF(m_sub_epair, SandrockSoedingreksoRhode)
    EPAIR_DEF(m_sub_epair, ForElectronPositron)

    py::class_<crosssection::EpairLPM, std::shared_ptr<crosssection::EpairLPM>>(
        m_sub_epair, "epair_lpm")
        .def(py::init<const ParticleDef&, const Medium&, double>(),
            py::arg("particle_def"), py::arg("medium"),
            py::arg("density_correction") = 1.0)
        .def("supression_factor", &crosssection::EpairLPM::suppression_factor,
            py::arg("energy"), py::arg("v"), py::arg("r2"), py::arg("beta"),
            py::arg("xi"));

    // --------------------------------------------------------------------- //
    // Annihilation
    // --------------------------------------------------------------------- //

    py::module m_sub_annihilation = m_sub.def_submodule("annihilation");
    py::class_<crosssection::Annihilation,
        std::shared_ptr<crosssection::Annihilation>,
        crosssection::Parametrization<Component>>(m_sub_annihilation, "Annihilation",
        R"pbdoc(

            Virtual class for the annihilation parametrizations. They can be initialized by using one of the given parametrizations with the following parameters

            Args:
                particle_def (:meth:`~proposal.particle.ParticleDef`): includes all static particle information for the parametrization such as mass, charge, etc.
                medium (:meth:`~proposal.medium`): includes all medium information for the parametrization such as densities or nucleon charges
                multiplier (double): Use a multiplicative factor for the differential crosssection. Can be used for testing or other studies

            The following parametrizations are currently implemented:

            * AnnihilationHeitler

            Example:
                To create a annihilation parametrization

                >>> positron = proposal.particle.EPlusDef.get()
                >>> medium = proposal.medium.StandardRock(1.0)
                >>> param = proposal.parametrization.annihilation.Heitler(positron, medium, 1.0)
                )pbdoc");

    py::class_<crosssection::AnnihilationHeitler,
        std::shared_ptr<crosssection::AnnihilationHeitler>,
        crosssection::Annihilation>(m_sub_annihilation, "Heitler")
        .def(py::init<>());

    // ---------------------------------------------------------------------
    // // Mupair
    // ---------------------------------------------------------------------
    // //

    py::module m_sub_mupair = m_sub.def_submodule("mupairproduction");
    py::class_<crosssection::MupairProduction,
        std::shared_ptr<crosssection::MupairProduction>,
        crosssection::Parametrization<Component>>(m_sub_mupair, "MupairProduction",
        R"pbdoc(

            Virtual class for the muon pair production parametrizations. They can be initialized by using one of the given parametrizations with the following parameters

            Args:
                particle_def (:meth:`~proposal.particle.ParticleDef`): includes all static particle information for the parametrization such as mass, charge, etc.
                medium (:meth:`~proposal.medium`): includes all medium information for the parametrization such as densities or nucleon charges
                energy_cuts (:meth:`~proposal.EnergyCutSettings`): energy cut setting for the parametrization
                multiplier (double): Use a multiplicative factor for the differential crosssection. Can be used for testing or other studies
                interpolation_def (:meth:`~proposal.InterpolationDef`): Only needed by Interpolant parametrizations. Includes settings for the interpolation
                particle_output (bool): If enabled, produced muons are sampled and returned as particles. Otherwise, only a DymamicData object is returned

            Since the differential cross section is given in :math:`\rho` as well, an intergration over this parameter is needed.
            When using the interpolation_def parameter, this integration is saved in interpolation tables (improving the performance of the calculation with neglible decline in accuracy).

            The following parametrizations are currently implemented:

            * KelnerKokoulinPetrukhin

            * KelnerKokoulinPetrukhinInterpolant

            Example:
                To create a muon pair production parametrization

                >>> mu = proposal.particle.MuMinusDef.get()
                >>> medium = proposal.medium.StandardRock(1.0)
                >>> cuts = proposal.EnergyCutSettings(-1, -1)
                >>> param = proposal.parametrization.mupairproduction.KelnerKokoulinPetrukhin(mu, medium, cuts, 1.0, True)
                )pbdoc");

    py::class_<crosssection::MupairProductionRhoIntegral,
        std::shared_ptr<crosssection::MupairProductionRhoIntegral>,
        crosssection::MupairProduction>(
        m_sub_mupair, "MupairProductionRhoIntegral")
        .def("function_to_integral",
            &crosssection::MupairProductionRhoIntegral::FunctionToIntegral,
             py::arg("particle_def"), py::arg("component"), py::arg("energy"),
             py::arg("v"), py::arg("rho"));

    MUPAIR_DEF(m_sub_mupair, KelnerKokoulinPetrukhin)

    // ---------------------------------------------------------------------
    // // Weak Interaction
    // ---------------------------------------------------------------------
    // //

    py::module m_sub_weak = m_sub.def_submodule("weakinteraction");
    py::class_<crosssection::WeakInteraction,
        std::shared_ptr<crosssection::WeakInteraction>,
        crosssection::Parametrization<Component>>(m_sub_weak, "WeakInteraction",
        R"pbdoc(

            Virtual class for the weak interaction parametrizations. They can be initialized by using one of the given parametrizations with the following parameters

            Args:
                particle_def (:meth:`~proposal.particle.ParticleDef`): includes all static particle information for the parametrization such as mass, charge, etc.
                medium (:meth:`~proposal.medium`): includes all medium information for the parametrization such as densities or nucleon charges
                multiplier (double): Use a multiplicative factor for the differential crosssection. Can be used for testing or other studies

            The following parametrizations are currently implemented:

            * WeakCooperSarkarMertsch

            Example:
                To create a weak interaction parametrization

                >>> mu = proposal.particle.MuMinusDef.get()
                >>> medium = proposal.medium.StandardRock(1.0)
                >>> param = proposal.parametrization.weakinteraction.WeakCooperSarkarMertsch(mu, medium, 1.0)
                )pbdoc");

    py::class_<crosssection::WeakCooperSarkarMertsch,
        std::shared_ptr<crosssection::WeakCooperSarkarMertsch>,
        crosssection::WeakInteraction>(m_sub_weak, "CooperSarkarMertsch")
        .def(py::init<>());

    // ---------------------------------------------------------------------
    // // Photo
    // ---------------------------------------------------------------------
    // //

    py::module m_sub_photo = m_sub.def_submodule("photonuclear");
    py::class_<crosssection::Photonuclear,
        std::shared_ptr<crosssection::Photonuclear>,
        crosssection::Parametrization<Component>>(m_sub_photo, "Photonuclear");

    // Shadow Effect
    py::class_<crosssection::ShadowEffect,
        std::shared_ptr<crosssection::ShadowEffect>>(m_sub_photo,
        "ShadowEffect",
        R"pbdoc(

            Virtual class for the parametrizations of the ShadowEffect used in the photonuclear interaction calculations.
            The nucleon shadowing describes the difference in the crosssection between the interaction of a photon with the whole nucleon compared to the interaction of a photon
            with a single nucleon (multplied by the number of nucleons in the atom). In general, the latter cross section in bigger than the real, measured crosssection, therefore
            this effect is called shadowing.

            The following parametrizations are currently implemented:

            * ShadowButkevichMikheyev

            * ShadowDuttaRenoSarcevicSeckel

                )pbdoc")
        .def("calculate_shadow_effect",
            &crosssection::ShadowEffect::CalculateShadowEffect,
            R"pbdoc(

            Calculate the shadow effect independently

            Args:
                Component (:meth:`~proposal.component`): Component to calculate the shadow effect
                x (float): Bjorken x
                nu (float): Fraction of energy transfered from the particle via the photon

                )pbdoc");

    py::class_<crosssection::ShadowDuttaRenoSarcevicSeckel,
        std::shared_ptr<crosssection::ShadowDuttaRenoSarcevicSeckel>,
        crosssection::ShadowEffect>(
        m_sub_photo, "ShadowDuttaRenoSarcevicSeckel")
        .def(py::init<>());

    py::class_<crosssection::ShadowButkevichMikheyev,
        std::shared_ptr<crosssection::ShadowButkevichMikheyev>,
        crosssection::ShadowEffect>(m_sub_photo, "ShadowButkevichMikheyev")
        .def(py::init<>());

    // Real Photon
    py::class_<crosssection::RealPhoton,
        std::shared_ptr<crosssection::RealPhoton>>(m_sub_photo, "RealPhoton")
        .def("calculate_hard_component",
            &crosssection::RealPhoton::CalculateHardComponent);

    py::class_<crosssection::SoftComponent,
        std::shared_ptr<crosssection::SoftComponent>, crosssection::RealPhoton>(
        m_sub_photo, "SoftComponent")
        .def(py::init<>());
    py::class_<crosssection::HardComponent,
        std::shared_ptr<crosssection::HardComponent>>(
        m_sub_photo, "HardComponent")
        .def(py::init<const ParticleDef&>(), py::arg("particle_def"))
        .def("calculate_hard_component", &crosssection::HardComponent::CalculateHardComponent);

    py::class_<crosssection::PhotoRealPhotonAssumption,
        std::shared_ptr<crosssection::PhotoRealPhotonAssumption>,
        crosssection::Photonuclear>(m_sub_photo, "PhotoRealPhotonAssumption",
        R"pbdoc(

            Virtual class for the parametrizations of photonuclear interaction. They can be initialized by using one of the given parametrizations with the following parameters

            Args:
                particle_def (:meth:`~proposal.particle.ParticleDef`): includes all static particle information for the parametrization such as mass, charge, etc.
                medium (:meth:`~proposal.medium`): includes all medium information for the parametrization such as densities or nucleon charges
                energy_cuts (:meth:`~proposal.EnergyCutSettings`): energy cut setting for the parametrization
                multiplier (double): Use a multiplicative factor for the differential crosssection. Can be used for testing or other studies
                hard_component (bool): Enabling or disabling the calculation of the hard component

            The following parametrizations are currently implemented:

            * Zeus

            * BezrukovBugaev

            * Kokoulin

            * Rhode

            Example:
                To create a photonuclear parametrization

                >>> mu = proposal.particle.MuMinusDef.get()
                >>> medium = proposal.medium.StandardRock(1.0)
                >>> cuts = proposal.EnergyCutSettings(-1, -1)
                >>> param = proposal.parametrization.photonuclear.Rhode(mu, medium, cuts, 1.0, True)
                )pbdoc");

    py::class_<crosssection::PhotoQ2Integral,
        std::shared_ptr<crosssection::PhotoQ2Integral>,
        crosssection::Photonuclear>(m_sub_photo, "PhotoQ2Integral",
        R"pbdoc(

            Virtual class for the parametrizations of photonuclear interaction. They can be initialized by using one of the given parametrizations with the following parameters

            Args:
                particle_def (:meth:`~proposal.particle.ParticleDef`): includes all static particle information for the parametrization such as mass, charge, etc.
                medium (:meth:`~proposal.medium`): includes all medium information for the parametrization such as densities or nucleon charges
                energy_cuts (:meth:`~proposal.EnergyCutSettings`): energy cut setting for the parametrization
                multiplier (double): Use a multiplicative factor for the differential crosssection. Can be used for testing or other studies
                ShadowEffect (:meth:`~proposal.parametrization.photonuclear.ShadowEffect`): Parametrization of the ShadowEffect to be used
                InterpolationDef (:meth:`~proposal.InterpolationDef`): Only needed by Interpolant parametrizations. Includes settings for the interpolation



            This virtual class includes all photonuclear interactions where the differential crosssection in given in :math:`Q^2`.
            The following parametrizations are currently implemented:

            * AbramowiczLevinLevyMaor91

            * AbramowiczLevinLevyMaor97

            * ButkevichMikheyev

            * RenoSarcevicSu

            * AbtFT

            * BlockDurandHa

            * AbramowiczLevinLevyMaor91Interpolant

            * AbramowiczLevinLevyMaor97Interpolant

            * ButkevichMikheyevInterpolant

            * RenoSarcevicSuInterpolant

            * AbtFTInterpolant

            * BlockDurandHaInterpolant

            The parametrization with "Interpolant" as a suffix creates an interpolation table for the :math:`Q^2` integration, which improved the perfomance.

            Example:
                To create a photonuclear parametrization

                >>> mu = proposal.particle.MuMinusDef.get()
                >>> medium = proposal.medium.StandardRock(1.0)
                >>> cuts = proposal.EnergyCutSettings(-1, -1)
                >>> shadow = proposal.parametrization.photonuclear.ShadowButkevichMikheyev()
                >>> interpol = proposal.InterpolationDef
                >>> param = proposal.parametrization.photonuclear.RenoSarcevicSuInterpolant(mu, medium, cuts, 1.0, shadow, interpol)
                )pbdoc");

    PHOTO_REAL_DEF(m_sub_photo, Zeus, RealPhotonAssumption)
    PHOTO_REAL_DEF(m_sub_photo, BezrukovBugaev, RealPhotonAssumption)
    PHOTO_REAL_DEF(m_sub_photo, Rhode, RealPhotonAssumption)
    PHOTO_REAL_DEF(m_sub_photo, Kokoulin,
        BezrukovBugaev) // Kokoulin derives from BezrukovBugaev

    PHOTO_Q2_DEF(m_sub_photo, AbramowiczLevinLevyMaor91)
    PHOTO_Q2_DEF(m_sub_photo, AbramowiczLevinLevyMaor97)
    PHOTO_Q2_DEF(m_sub_photo, ButkevichMikheyev)
    PHOTO_Q2_DEF(m_sub_photo, RenoSarcevicSu)
    PHOTO_Q2_DEF(m_sub_photo, AbtFT)
    PHOTO_Q2_DEF(m_sub_photo, BlockDurandHa)

    // --------------------------------------------------------------------- //
    // Ionization
    // --------------------------------------------------------------------- //

    py::module m_sub_ioniz = m_sub.def_submodule("ionization");

    py::class_<crosssection::Ionization,
        std::shared_ptr<crosssection::Ionization>,
        crosssection::Parametrization<Medium>>(m_sub_ioniz, "Ionization",
        R"pbdoc(

            Virtual class for the Ionization parametrizations. They can be initialized by using one of the given parametrizations with the following parameters

            Args:
                particle_def (:meth:`~proposal.particle.ParticleDef`): includes all static particle information for the parametrization such as mass, charge, etc.
                medium (:meth:`~proposal.medium`): includes all medium information for the parametrization such as densities or nucleon charges
                energy_cuts (:meth:`~proposal.EnergyCutSettings`): energy cut setting for the parametrization
                multiplier (double): Use a multiplicative factor for the differential crosssection. Can be used for testing or other studies

            The following parametrizations are currently implemented:

            * BetheBlochRossi

            * IonizBergerSeltzerBhabha (for positron propagation)

            * IonizBergerSeltzerMoller (for electron propagation)

            Example:
                To create a ionization parametrization

                >>> mu = proposal.particle.MuMinusDef.get()
                >>> medium = proposal.medium.StandardRock(1.0)
                >>> cuts = proposal.EnergyCutSettings(-1, -1)
                >>> param = proposal.parametrization.ionization.BetheBlochRossi(mu, medium, cuts, multiplier)
                )pbdoc")
        .def("differential_crosssection",
            py::overload_cast<ParticleDef const&, Medium const&, double,
                double>(&crosssection::Ionization::DifferentialCrossSection,
                py::const_),
            py::arg("particle_def"), py::arg("medium"), py::arg("energy"),
            py::arg("v"))
        .def("kinematic_limits",
            py::overload_cast<ParticleDef const&, Medium const&, double>(
                &crosssection::Ionization::GetKinematicLimits, py::const_),
            py::arg("particle_def"), py::arg("medium"), py::arg("energy"));

    py::class_<crosssection::IonizBetheBlochRossi,
        std::shared_ptr<crosssection::IonizBetheBlochRossi>,
        crosssection::Ionization>(m_sub_ioniz, "BetheBlochRossi")
        .def(py::init<const EnergyCutSettings&>(), py::arg("energy_cuts"));

    py::class_<crosssection::IonizBergerSeltzerBhabha,
        std::shared_ptr<crosssection::IonizBergerSeltzerBhabha>,
        crosssection::Ionization>(m_sub_ioniz, "BergerSeltzerBhabha")
        .def(py::init<const EnergyCutSettings&>(), py::arg("energy_cuts"));

    py::class_<crosssection::IonizBergerSeltzerMoller,
        std::shared_ptr<crosssection::IonizBergerSeltzerMoller>,
        crosssection::Ionization>(m_sub_ioniz, "BergerSeltzerMoller")
        .def(py::init<const EnergyCutSettings&>(), py::arg("energy_cuts"));

    // Photon interactions

    // --------------------------------------------------------------------- //
    // Compton Scattering
    // --------------------------------------------------------------------- //

    py::module m_sub_compton = m_sub.def_submodule("compton");

    py::class_<crosssection::Compton, std::shared_ptr<crosssection::Compton>,
        crosssection::Parametrization<Component>>(m_sub_compton, "Compton",
        R"pbdoc(

            Virtual class for the Compton scattering parametrizations. They can be initialized by using one of the given parametrizations with the following parameters

            Args:
                particle_def (:meth:`~proposal.particle.ParticleDef`): includes all static particle information for the parametrization such as mass, charge, etc.
                medium (:meth:`~proposal.medium`): includes all medium information for the parametrization such as densities or nucleon charges
                energy_cuts (:meth:`~proposal.EnergyCutSettings`): energy cut setting for the parametrization
                multiplier (double): Use a multiplicative factor for the differential crosssection. Can be used for testing or other studies

            The following parametrizations are currently implemented:

            * KleinNishina

            Example:
                To create a compton scattering parametrization

                >>> gamma = proposal.particle.GammaDef.get()
                >>> medium = proposal.medium.StandardRock(1.0)
                >>> cuts = proposal.EnergyCutSettings(-1, -1)
                >>> param = proposal.parametrization.compton.KleinNishina(gamma, medium, cuts, multiplier)
                )pbdoc");

    py::class_<crosssection::ComptonKleinNishina,
        std::shared_ptr<crosssection::ComptonKleinNishina>,
        crosssection::Compton>(m_sub_compton, "KleinNishina")
        .def(py::init<>());

    // --------------------------------------------------------------------- //
    // PhotoPairProduction
    // --------------------------------------------------------------------- //

    py::module m_sub_photopair = m_sub.def_submodule("photopair");

    py::class_<crosssection::PhotoPairProduction,
        std::shared_ptr<crosssection::PhotoPairProduction>,
        crosssection::Parametrization<Component>>(m_sub_photopair, "PhotoPair",
        R"pbdoc(

            Virtual class for the PhotoPairProduction parametrizations. They can be initialized by using one of the given parametrizations with the following parameters

            Args:
                particle_def (:meth:`~proposal.particle.ParticleDef`): includes all static particle information for the parametrization such as mass, charge, etc.
                medium (:meth:`~proposal.medium`): includes all medium information for the parametrization such as densities or nucleon charges
                multiplier (double): Use a multiplicative factor for the differential crosssection. Can be used for testing or other studies

            The following parametrizations are currently implemented:

            * Tsai

            Example:
                To create a photopair parametrization

                >>> gamma = proposal.particle.GammaDef.get()
                >>> medium = proposal.medium.StandardRock(1.0)
                >>> param = proposal.parametrization.PhotoPair.Tsai(gamma, medium, multiplier)
                )pbdoc");

    py::class_<crosssection::PhotoPairTsai,
        std::shared_ptr<crosssection::PhotoPairTsai>,
        crosssection::PhotoPairProduction>(m_sub_photopair, "Tsai")
        .def(py::init<>());

    py::class_<crosssection::KinematicLimits,
    std::shared_ptr<crosssection::KinematicLimits>>(m_sub, "KinematicLimits")
        .def(py::init<>())
        .def_readwrite("v_min", &crosssection::KinematicLimits::v_min)
        .def_readwrite("v_max", &crosssection::KinematicLimits::v_max)
        .def("__repr__",
         [](const crosssection::KinematicLimits &lim) {
             return "(v_min: " + std::to_string(lim.v_min) + ", v_max: " + std::to_string(lim.v_max) + ")";
         });
    // py::class_<PhotoAngleDistribution,
    // std::shared_ptr<PhotoAngleDistribution>>(
    //     m_sub_photopair, "PhotoAngleDistribution",
    //     R"pbdoc(

    //         Virtual class for the PhotoAngleDistribution parametrizations.
    //         They can be initialized by using one of the given
    //         parametrizations with the following parameters

    //         Args:
    //             particle_def (:meth:`~proposal.particle.ParticleDef`):
    //             includes all static particle information for the
    //             parametrization such as mass, charge, etc. medium
    //             (:meth:`~proposal.medium`): includes all medium information
    //             for the parametrization such as densities or nucleon charges

    //         The following parametrizations are currently implemented:

    //         * PhotoPairNoDeflection

    //         * PhotoPairTsaiIntegral

    //         * PhotoPairEGS

    //         Example:
    //             To create a PhotoAngleDistribution parametrization

    //             >>> gamma = proposal.particle.GammaDef.get()
    //             >>> medium = proposal.medium.StandardRock(1.0)
    //             >>> param =
    //             proposal.parametrization.PhotoAngleDistribution.TsaiIntegral(gamma,
    //             medium) )pbdoc")
    //     .def("SetCurrentComponent",
    //         &PhotoAngleDistribution::SetCurrentComponent,
    //         py::arg("component_index"))
    //     .def("SampleAngles", &PhotoAngleDistribution::SampleAngles,
    //         py::arg("energy"), py::arg("rho"), py::arg("component_index"));

    // py::class_<PhotoAngleTsaiIntegral,
    // std::shared_ptr<PhotoAngleTsaiIntegral>,
    //     PhotoAngleDistribution>(m_sub_photopair, "PhotoAngleTsaiIntegral")
    //     .def(py::init<const ParticleDef&, std::shared_ptr<const Medium>>(),
    //         py::arg("particle_def"), py::arg("medium"))
    //     .def("FunctionToIntegral",
    //     &PhotoAngleTsaiIntegral::FunctionToIntegral,
    //         py::arg("energy"), py::arg("x"), py::arg("theta"));

    // py::class_<PhotoAngleNoDeflection,
    // std::shared_ptr<PhotoAngleNoDeflection>,
    //     PhotoAngleDistribution>(m_sub_photopair, "PhotoAngleNoDeflection")
    //     .def(py::init<const ParticleDef&, std::shared_ptr<const Medium>>(),
    //         py::arg("particle_def"), py::arg("medium"));

    // py::class_<PhotoAngleEGS, std::shared_ptr<PhotoAngleEGS>,
    //     PhotoAngleDistribution>(m_sub_photopair, "PhotoAngleEGS")
    //     .def(py::init<const ParticleDef&, std::shared_ptr<const Medium>>(),
    //         py::arg("particle_def"), py::arg("medium"));

    // py::class_<PhotoAngleDistribution::DeflectionAngles,
    //     std::shared_ptr<PhotoAngleDistribution::DeflectionAngles>>(
    //     m_sub_photopair, "DeflectionAngles")
    //     .def(py::init<>())
    //     .def_readwrite(
    //         "cosphi0", &PhotoAngleDistribution::DeflectionAngles::cosphi0)
    //     .def_readwrite(
    //         "theta0", &PhotoAngleDistribution::DeflectionAngles::theta0)
    //     .def_readwrite(
    //         "cosphi1", &PhotoAngleDistribution::DeflectionAngles::cosphi1)
    //     .def_readwrite(
    //         "theta1", &PhotoAngleDistribution::DeflectionAngles::theta1);
}
