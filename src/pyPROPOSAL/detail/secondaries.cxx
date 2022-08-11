#include "PROPOSAL/secondaries/parametrization/Parametrization.h"
#include "PROPOSAL/secondaries/parametrization/annihilation/HeitlerAnnihilation.h"
#include "PROPOSAL/secondaries/parametrization/bremsstrahlung/BremsNoDeflection.h"
#include "PROPOSAL/secondaries/parametrization/bremsstrahlung/BremsEGS4Approximation.h"
#include "PROPOSAL/secondaries/parametrization/bremsstrahlung/BremsKochMotz.h"
#include "PROPOSAL/secondaries/parametrization/compton/NaivCompton.h"
#include "PROPOSAL/secondaries/parametrization/epairproduction/KelnerKokoulinPetrukhinEpairProduction.h"
#include "PROPOSAL/secondaries/parametrization/epairproduction/NaivEpairProduction.h"
#include "PROPOSAL/secondaries/parametrization/ionization/NaivIonization.h"
#include "PROPOSAL/secondaries/parametrization/photomupairproduction/PhotoMuPairProductionBurkhardtKelnerKokoulin.h"
#include "PROPOSAL/secondaries/parametrization/mupairproduction/KelnerKokoulinPetrukhinMupairProduction.h"
#include "PROPOSAL/secondaries/parametrization/photonuclear/Photonuclear.h"
#include "PROPOSAL/secondaries/parametrization/photoproduction/Photoproduction.h"
#include "PROPOSAL/secondaries/parametrization/photoeffect/Photoeffect.h"
#include "PROPOSAL/secondaries/parametrization/photoeffect/PhotoeffectNoDeflection.h"
#include "PROPOSAL/secondaries/parametrization/photopairproduction/PhotoPairProductionTsai.h"
#include "PROPOSAL/secondaries/parametrization/photopairproduction/PhotoPairProductionKochMotz.h"
#include "PROPOSAL/secondaries/parametrization/weakinteraction/WeakCooperSarkarMertsch.h"

#include "PROPOSAL/secondaries/SecondariesCalculator.h"

#include "pyPROPOSAL/pyBindings.h"

namespace py = pybind11;
using namespace PROPOSAL;

template <typename T, typename BaseT> struct SecondariesBuilder {
    template <typename... Args> auto decl_virtual_param(Args... args)
    {
        py::class_<T, BaseT, std::shared_ptr<T>>(args...);
    };

    template <typename... Args> auto decl_param(Args... args)
    {
        py::class_<T, BaseT, std::shared_ptr<T>>(args...).def(
            py::init<ParticleDef, Medium>(), py::arg("particle_def"),
            py::arg("medium"));
    }

    template <typename... Args> auto decl_rho_param(Args... args)
    {
        py::class_<T, BaseT, std::shared_ptr<T>>(args...)
            .def(py::init<ParticleDef, Medium>(), py::arg("particle_def"),
                py::arg("medium"))
            .def("calculate_rho", &T::CalculateRho,
                 R"pbdoc(
                Samples the asymmetry factor for interactions where two particles
                are created. For EpairProduction and MupairProduction, this
                factor is defined as :math:`\frac{E_+ - E_-}{E_+ + E_-}`, where
                :math:`E_-` is the energy of the created particle and :math:`E_+`
                the energy of the created antiparticle. For annihilation and
                photopairproduction, this factor is defined as
                :math:`\frac{E_{\gamma,1}}{E_+ + m_e}`, respectively
                :math:`\frac{E_-}{E_{\gamma}}`.
            )pbdoc");
    }
};

void init_secondaries(py::module& m)
{
    py::module m_sub = m.def_submodule("secondaries");

    py::class_<secondaries::Parametrization,
        std::shared_ptr<secondaries::Parametrization>>(m_sub, "Parametrization")
        .def_property_readonly("random_numbers",
            &secondaries::Parametrization::RequiredRandomNumbers)
        .def("interaction_type",
            &secondaries::Parametrization::GetInteractionType)
        .def("calculate_secondaries",
            &secondaries::Parametrization::CalculateSecondaries);

    SecondariesBuilder<secondaries::Annihilation, secondaries::Parametrization> {}
        .decl_virtual_param(m_sub, "Annihilation");
    SecondariesBuilder<secondaries::Bremsstrahlung,
        secondaries::Parametrization> {}
        .decl_virtual_param(m_sub, "Bremsstrahlung");
    SecondariesBuilder<secondaries::Compton, secondaries::Parametrization> {}
        .decl_virtual_param(m_sub, "Compton");
    SecondariesBuilder<secondaries::EpairProduction,
        secondaries::Parametrization> {}
        .decl_virtual_param(m_sub, "EpairProduction");
    SecondariesBuilder<secondaries::Ionization, secondaries::Parametrization> {}
        .decl_virtual_param(m_sub, "Ionization");
    SecondariesBuilder<secondaries::MupairProduction,
        secondaries::Parametrization> {}
        .decl_virtual_param(m_sub, "MupairProduction");
    SecondariesBuilder<secondaries::PhotoMuPairProduction,
        secondaries::Parametrization> {}
        .decl_virtual_param(m_sub, "PhotoMuPairProduction");
    SecondariesBuilder<secondaries::Photonuclear,
        secondaries::Parametrization> {}
        .decl_param(m_sub, "Photonuclear");
    SecondariesBuilder<secondaries::Photoproduction,
        secondaries::Parametrization> {}
        .decl_param(m_sub, "Photoproduction");
    SecondariesBuilder<secondaries::Photoeffect,
        secondaries::Parametrization> {}
        .decl_virtual_param(m_sub, "Photoeffect");
    SecondariesBuilder<secondaries::PhotoPairProduction,
        secondaries::Parametrization> {}
        .decl_virtual_param(m_sub, "PhotoPairProduction");
    SecondariesBuilder<secondaries::WeakInteraction,
        secondaries::Parametrization> {}
        .decl_virtual_param(m_sub, "WeakInteraction");

    SecondariesBuilder<secondaries::HeitlerAnnihilation,
        secondaries::Annihilation> {}
        .decl_rho_param(m_sub, "HeitlerAnnihilation");

    SecondariesBuilder<secondaries::BremsNoDeflection,
        secondaries::Bremsstrahlung> {}
        .decl_param(m_sub, "BremsNoDeflection");
    SecondariesBuilder<secondaries::BremsEGS4Approximation,
        secondaries::Bremsstrahlung> {}
        .decl_param(m_sub, "BremsEGS4Approximation");
    SecondariesBuilder<secondaries::BremsKochMotz,
        secondaries::Bremsstrahlung> {}
        .decl_param(m_sub, "BremsKochMotz");

    SecondariesBuilder<secondaries::NaivCompton, secondaries::Compton> {}
        .decl_param(m_sub, "NaivCompton");

    SecondariesBuilder<secondaries::KelnerKokoulinPetrukhinEpairProduction,
        secondaries::EpairProduction> {}
        .decl_rho_param(m_sub, "KelnerKokoulinPetrukhinEpairProduction");
    SecondariesBuilder<secondaries::NaivEpairProduction,
        secondaries::EpairProduction> {}
        .decl_param(m_sub, "NaivEpairProduction");

    SecondariesBuilder<secondaries::NaivIonization, secondaries::Ionization> {}
        .decl_param(m_sub, "NaivIonization");

    SecondariesBuilder<secondaries::KelnerKokoulinPetrukhinMupairProduction,
        secondaries::MupairProduction> {}
        .decl_rho_param(m_sub, "KelnerKokoulinPetrukhinMupairProduction");

    py::class_<secondaries::PhotoMuPairProductionBurkhardtKelnerKokoulin, secondaries::PhotoMuPairProduction,
    std::shared_ptr<secondaries::PhotoMuPairProductionBurkhardtKelnerKokoulin>>(m_sub,
            "PhotoMuPairProductionBurkhardtKelnerKokoulin")
        .def(py::init<ParticleDef, Medium>(), py::arg("particle_def"),
             py::arg("medium"))
        .def("calculate_x", &secondaries::PhotoMuPairProductionBurkhardtKelnerKokoulin::Calculatex,
             py::arg("energy"), py::arg("rnd"), py::arg("component"));

    SecondariesBuilder<secondaries::PhotoPairProductionTsai,
        secondaries::PhotoPairProduction> {}
        .decl_rho_param(m_sub, "PhotoTsai");
    SecondariesBuilder<secondaries::PhotoPairProductionTsaiForwardPeaked,
        secondaries::PhotoPairProduction> {}
        .decl_rho_param(m_sub, "PhotoTsaiForwardPeaked");
    SecondariesBuilder<secondaries::PhotoPairProductionKochMotzForwardPeaked,
        secondaries::PhotoPairProduction> {}
        .decl_rho_param(m_sub, "PhotoKochMotzForwardPeaked");

    SecondariesBuilder<secondaries::PhotoeffectNoDeflection, secondaries::Photoeffect> {}
        .decl_param(m_sub, "PhotoeffectNoDeflection");

    py::class_<secondaries::WeakCooperSarkarMertsch, secondaries::WeakInteraction,
        std::shared_ptr<secondaries::WeakCooperSarkarMertsch>>(m_sub,
                "WeakCooperSarkarMertsch")
        .def(py::init<ParticleDef, Medium>(), py::arg("particle_def"),
             py::arg("medium"))
        .def("calculate_relative_loss", &secondaries::WeakCooperSarkarMertsch::CalculateRelativeLoss,
             py::arg("energy"), py::arg("rnd"), py::arg("component"));

    py::class_<SecondariesCalculator, std::shared_ptr<SecondariesCalculator>>(
        m_sub, "SecondariesCalculator")
        .def(py::init<std::vector<InteractionType>, ParticleDef const&,
            Medium const&>())
        .def("n_rnd", &SecondariesCalculator::RequiredRandomNumbers)
        .def("calculate_secondaries",
            &SecondariesCalculator::CalculateSecondaries);

    m_sub.def("make_secondary",
        [](InteractionType type, ParticleDef const& p, Medium const& m) {
            return std::shared_ptr<secondaries::Parametrization>(
                DefaultFactory<secondaries::Parametrization>::Create(
                    type, p, m));
        });
}
