#include "PROPOSAL/PROPOSAL.h"
#include "pyBindings.h"

namespace py = pybind11;
using namespace PROPOSAL;

void init_crosssection(py::module& m) {
    py::module m_sub = m.def_submodule("crosssection");

    py::class_<CrossSection, std::shared_ptr<CrossSection>>(m_sub,
                                                            "CrossSection",
                                                            R"pbdoc( 

            Virtual class for crosssections. The crosssection class provides all mathematical methods to process the theoretical, differential crosssections that are
            given by the parametrizations. A cross section class can be initialized with the following parameters

            Args:                                                                                                  
                param (:meth:`~proposal.parametrization`): parametrization for the crosssection, including the chosen theoretical model
                interpolation_def (:meth:`~proposal.InterpolationDef`): Only needed by Interpolant parametrizations. Includes settings for the interpolation                                    
 
            The crosssection class can either work with interpolation tables or with exact intergration for every single calculation.
            Since the usage of interpolation tables can improve the speed of the propagation by several orders of magnitude (with neglible decline in accuracy) it is highly recommended
            to use the interpolation methods.
                                                            
            There are specific crosssection classes for every interaction that can be used. 

            For propagation of massive leptons, there are the following crosssections:

            * AnnihilationIntegral / AnnihilationInterpolant

            * BremsIntegral / BremsInterpolant

            * EpairIntegral / EpairInterpolant

            * IonizIntegral / IonizInterpolant

            * MupairIntegral / MupairInterpolant

            * PhotoIntegral / PhotoInterpolant

            * WeakIntegral / WeakInterpolant

            For propagation of photons, there are the following crosssections:

            * ComptonIntegral / ComptonInterpolant

            * PhotoPairIntegral / PhotoPairInterpolant

            Example:
                To create a bremsstrahlung CrossSection

                >>> mu = proposal.particle.MuMinusDef.get()
                >>> medium = proposal.medium.StandardRock(1.0)
                >>> cuts = proposal.EnergyCutSettings(-1, -1)
                >>> interpol = proposal.InterpolationDef
                >>> param = proposal.parametrization.bremsstrahlung.SandrockSoedingreksoRhode(mu, medium, cuts, 1.0, False)
                >>> cross = proposal.crosssection.BremsInterpolant(param, interpol)
                >>> cross.calculate_dEdx(1e6) # exmaple usage of the created crosssection class...
                )pbdoc")
        .def("__str__", &py_print<CrossSection>)
        .def("calculate_dEdx", &CrossSection::CalculatedEdx, py::arg("energy"),
             R"pbdoc( 

            Calculates the continous energy loss :math:`\langle \frac{dE}{dx} \rangle`, which equals to

                .. math:: \frac{N_A}{A} \cdot E \cdot \int_{v_{min}}^{v_{cut}} v \cdot \frac{d\sigma}{dv} dv

            with the particle energy E, the relative energy loss v and the crosssection :math:`\sigma`. The value :math:`v_{cut}` is the energy cut to differentiate between
            continous and stochastic losses in PROPOSAL, see :meth:`~proposal.EnergyCutSettings` for more information on the energy cuts.

            Args:                                                                                                  
                energy (float): energy in MeV
                                                            
                )pbdoc")
        .def("calculate_dE2dx", &CrossSection::CalculatedE2dx,
             py::arg("energy"),
             R"pbdoc( 

            Calculates the value

                .. math:: \frac{N_A}{A} \cdot E^2 \cdot \int_{v_{min}}^{v_{cut}} v^2 \cdot \frac{d\sigma}{dv} dv

            with the particle energy E, the relative energy loss v and the crosssection :math:`\sigma`. The value :math:`v_{cut}` is the energy cut to differentiate between
            continous and stochastic losses in PROPOSAL, see :meth:`~proposal.EnergyCutSettings` for more information on the energy cuts.

            The value is important for the calculation of the ContinuousRandomization (see :meth:`~proposal.ContinuousRandomizer`) 

            Args:                                                                                                  
                energy (float): energy in MeV
                                                            
                )pbdoc")
        .def("calculate_dNdx",
             (double (CrossSection::*)(double)) & CrossSection::CalculatedNdx,
             py::arg("energy"),
             R"pbdoc( 

            Calculates the total cross section

                .. math:: \frac{N_A}{A} \cdot \int_{v_{cut}}^{v_{max}} \frac{d\sigma}{dv} dv

            with the particle energy E, the relative energy loss v and the crosssection :math:`\sigma`. The value v_{cut} is the energy cut to differentiate between
            continous and stochastic losses in PROPOSAL, see :meth:`~proposal.EnergyCutSettings` for more information on the energy cuts.

            Note that this integral only includes the v values about our cut, therefore this values represents only the total crosssection for the stochastic energy losses. 

            Args:                                                                                                  
                energy (float): energy in MeV
                                                            
                )pbdoc")
        .def("calculate_dNdx_rnd",
             (double (CrossSection::*)(double, double)) &
                 CrossSection::CalculatedNdx,
             py::arg("energy"), py::arg("rnd"),

             R"pbdoc( 

            Calculates the total cross section

                .. math:: \frac{N_A}{A} \cdot \int_{v_{cut}}^{v_{max}} \frac{d\sigma}{dv} dv

            with the particle energy E, the relative energy loss v and the crosssection :math:`\sigma`. The value v_{cut} is the energy cut to differentiate between
            continous and stochastic losses in PROPOSAL, see :meth:`~proposal.EnergyCutSettings` for more information on the energy cuts.

            Furthermore, for every component in the medium, a stochatic energy loss in saved based on the random number rnd. The values are saved in the crosssection class and can be used
            by other methods.

            Note that this integral only includes the v values about our cut, therefore this values represents only the total crosssection for the stochastic energy losses. 

            Args:                                                                                                  
                energy (float): energy in MeV
                rnd (float): random number between 0 and 1
                                               
                )pbdoc")
        .def("calculate_stochastic_loss",
             (double (CrossSection::*)(double, double, double)) &
                 CrossSection::CalculateStochasticLoss,
             py::arg("energy"), py::arg("rnd1"), py::arg("rnd2"),
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
        .def("calculate_produced_particles",
             &CrossSection::CalculateProducedParticles, py::arg("energy"),
             py::arg("energy_loss"), py::arg("initial_direction"),
             R"pbdoc( 

            If particles are produced in the interaction of this CrossSection, the methods samples those particles corresponding to the energy of the initial particle as
            well as the energy loss of the initial particle

            Args:                                                                                                  
                energy (float): primary particle energy in MeV
                energy_loss (float): energy loss of the primary particle in MeV
                initial_direction (Vector3D): direction of the parent particle

            Returns:
                List of created particles as well as a boolean with the information whether the initial particle has been destroyed in the interaction
                     
                )pbdoc")
        .def_property_readonly("id", &CrossSection::GetTypeId,
                               R"pbdoc( 

            Internal id of the current interaction, see :meth:`~proposal.DynamicData` for all available id's.
                                                            
                )pbdoc")
        .def_property_readonly("parametrization",
                               &CrossSection::GetParametrization,
                               R"pbdoc( 

            Pointer to the current parametrization object
                                                            
                )pbdoc")        
        .def("calculate_cumulative_crosssection", &CrossSection::CalculateCumulativeCrossSection,
             py::arg("energy"),
             py::arg("component"),
             py::arg("v"),
             R"pbdoc(

            Return the cumulative cross section, e.g.

                .. math:: \int_{v_{cut}}^{v} \frac{d\sigma}{dv} dv

            Args:
                energy (float): energy in MeV
                component (int): number of the component to evaluate the crosssection. Use component = 0 for media with only one component
                v (float): limit in v to evaluate the differential cross section

            Returns:
                cumulative crosssection, e.g. the probability for a stochastic loss where v is smaller than the given parameter (and bigger than the defined v_cut)
                )pbdoc");

    py::class_<CrossSectionIntegral, std::shared_ptr<CrossSectionIntegral>,
               CrossSection>(m_sub, "CrossSectionIntegral");
    py::class_<CrossSectionInterpolant,
               std::shared_ptr<CrossSectionInterpolant>, CrossSection>(
        m_sub, "CrossSectionInterpolant");

    py::class_<BremsIntegral, std::shared_ptr<BremsIntegral>, CrossSectionIntegral>(m_sub, "BremsIntegral")
        .def(py::init<const Bremsstrahlung&>(), py::arg("parametrization"));
    py::class_<EpairIntegral, std::shared_ptr<EpairIntegral>, CrossSectionIntegral>(m_sub, "EpairIntegral")
        .def(py::init<const EpairProduction&>(), py::arg("parametrization"));
    py::class_<PhotoIntegral, std::shared_ptr<PhotoIntegral>, CrossSectionIntegral>(m_sub, "PhotoIntegral")
        .def(py::init<const Photonuclear&>(), py::arg("parametrization"));
    py::class_<IonizIntegral, std::shared_ptr<IonizIntegral>, CrossSectionIntegral>(m_sub, "IonizIntegral")
        .def(py::init<const Ionization&>(), py::arg("parametrization"));
    py::class_<AnnihilationIntegral, std::shared_ptr<AnnihilationIntegral>, CrossSectionIntegral>(m_sub, "AnnihilationIntegral")
        .def(py::init<const Annihilation&>(), py::arg("parametrization"));
    py::class_<MupairIntegral, std::shared_ptr<MupairIntegral>, CrossSectionIntegral>(m_sub, "MupairIntegral")
        .def(py::init<const MupairProduction&>(), py::arg("parametrization"));
    py::class_<WeakIntegral, std::shared_ptr<WeakIntegral>, CrossSectionIntegral>(m_sub, "WeakIntegral")
         .def(py::init<const WeakInteraction&>(), py::arg("parametrization"));
    py::class_<ComptonIntegral, std::shared_ptr<ComptonIntegral>, CrossSectionIntegral>(m_sub, "ComptonIntegral")
         .def(py::init<const Compton&>(), py::arg("parametrization"));
    py::class_<PhotoPairIntegral, std::shared_ptr<PhotoPairIntegral>, CrossSectionIntegral>(m_sub, "PhotoPairIntegral")
         .def(py::init<const PhotoPairProduction&, const PhotoAngleDistribution&>(), py::arg("parametrization"), py::arg("photoangle_distribution"));

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
    py::class_<AnnihilationInterpolant, std::shared_ptr<AnnihilationInterpolant>, CrossSectionInterpolant>(
        m_sub, "AnnihilationInterpolant")
        .def(py::init<const Annihilation&, InterpolationDef>(), py::arg("parametrization"), py::arg("interpolation_def"));
    py::class_<MupairInterpolant, std::shared_ptr<MupairInterpolant>, CrossSectionInterpolant>(
        m_sub, "MupairInterpolant")
        .def(py::init<const MupairProduction&, InterpolationDef>(), py::arg("parametrization"), py::arg("interpolation_def"));
    py::class_<WeakInterpolant, std::shared_ptr<WeakInterpolant>, CrossSectionInterpolant>(
        m_sub, "WeakInterpolant")
        .def(py::init<const WeakInteraction&, InterpolationDef>(), py::arg("parametrization"), py::arg("interpolation_def"));
    py::class_<ComptonInterpolant, std::shared_ptr<ComptonInterpolant>, CrossSectionInterpolant>(
        m_sub, "ComptonInterpolant")
        .def(py::init<const Compton&, InterpolationDef>(), py::arg("parametrization"), py::arg("interpolation_def"));
    py::class_<PhotoPairInterpolant, std::shared_ptr<PhotoPairInterpolant>, CrossSectionInterpolant>(
        m_sub, "PhotoPairInterpolant")
        .def(py::init<const PhotoPairProduction&, const PhotoAngleDistribution&, InterpolationDef>(), py::arg("parametrization"), py::arg("photoangle_distribution"), py::arg("interpolation_def"));
}
