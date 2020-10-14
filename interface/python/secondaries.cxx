#include "PROPOSAL/secondaries/Parametrization.h"

#include "PROPOSAL/secondaries/annihilation/SingleDifferentialAnnihilation.h"
#include "PROPOSAL/secondaries/bremsstrahlung/NaivBremsstrahlung.h"
#include "PROPOSAL/secondaries/compton/NaivCompton.h"
#include "PROPOSAL/secondaries/epairproduction/KelnerKokoulinPetrukhinEpairProduction.h"
#include "PROPOSAL/secondaries/ionization/NaivIonization.h"
#include "PROPOSAL/secondaries/mupairproduction/KelnerKokoulinPetrukhinMupairProduction.h"
#include "PROPOSAL/secondaries/photonuclear/Photonuclear.h"
#include "PROPOSAL/secondaries/photopairproduction/PhotoTsai.h"
#include "PROPOSAL/secondaries/weakinteraction/NaivWeakInteraction.h"

#include "PROPOSAL/secondaries/SecondariesCalculator.h"

#include "pyBindings.h"

namespace py = pybind11;
using namespace PROPOSAL;

void init_secondaries(py::module& m)
{
    py::module m_sub = m.def_submodule("secondaries");

    py::class_<secondaries::Parametrization, std::shared_ptr<secondaries::Parametrization>>( m_sub, "Parametrization")
        .def_property_readonly("random_numbers", &secondaries::Parametrization::RequiredRandomNumbers)
        .def("interaction_type", &secondaries::Parametrization::GetInteractionType)
        .def("calculate_secondaries", &secondaries::Parametrization::CalculateSecondaries);

    py::class_<secondaries::Annihilation, secondaries::Parametrization, std::shared_ptr<secondaries::Annihilation>>(m_sub, "Annihilation");
    py::class_<secondaries::SingleDifferentialAnnihilation, secondaries::Annihilation, std::shared_ptr<secondaries::SingleDifferentialAnnihilation>>(m_sub, "SingleDifferentialAnnihilation");

    py::class_<secondaries::Bremsstrahlung, secondaries::Parametrization, std::shared_ptr<secondaries::Bremsstrahlung>>(m_sub, "Bremsstrahlung");
    py::class_<secondaries::NaivBremsstrahlung, secondaries::Bremsstrahlung, std::shared_ptr<secondaries::NaivBremsstrahlung>>(m_sub, "NaivBremsstrahlung");

    py::class_<secondaries::Compton, secondaries::Parametrization, std::shared_ptr<secondaries::Compton>>(m_sub, "Compton");
    py::class_<secondaries::NaivCompton, secondaries::Compton, std::shared_ptr<secondaries::NaivCompton>>(m_sub, "NaivCompton");

    py::class_<secondaries::EpairProduction, secondaries::Parametrization, std::shared_ptr<secondaries::EpairProduction>>(m_sub, "EpairProduction");
    py::class_<secondaries::KelnerKokoulinPetrukhinEpairProduction, secondaries::EpairProduction, std::shared_ptr<secondaries::KelnerKokoulinPetrukhinEpairProduction>>(m_sub, "KelnerKokoulinPetrukhinEpairProduction");

    py::class_<secondaries::Ionization, secondaries::Parametrization, std::shared_ptr<secondaries::Ionization>>(m_sub, "Ionization");
    py::class_<secondaries::NaivIonization, secondaries::Ionization, std::shared_ptr<secondaries::NaivIonization>>(m_sub, "NaivIonization");

    py::class_<secondaries::MupairProduction, secondaries::Parametrization, std::shared_ptr<secondaries::MupairProduction>>(m_sub, "MupairProduction");
    py::class_<secondaries::KelnerKokoulinPetrukhinMupairProduction, secondaries::MupairProduction, std::shared_ptr<secondaries::KelnerKokoulinPetrukhinMupairProduction>>(m_sub, "KelnerKokoulinPetrukhinMupairProduction");

    py::class_<secondaries::Photonuclear, secondaries::Parametrization, std::shared_ptr<secondaries::Photonuclear>>(m_sub, "Photonuclear");

    py::class_<secondaries::PhotopairProduction, secondaries::Parametrization, std::shared_ptr<secondaries::PhotopairProduction>>(m_sub, "PhotopairProduction");
    py::class_<secondaries::PhotoTsai, secondaries::PhotopairProduction, std::shared_ptr<secondaries::PhotoTsai>>(m_sub, "PhotoTsai");

    py::class_<secondaries::WeakInteraction, secondaries::Parametrization, std::shared_ptr<secondaries::WeakInteraction>>(m_sub, "WeakInteraction");
    py::class_<secondaries::NaivWeakInteraction, secondaries::WeakInteraction, std::shared_ptr<secondaries::NaivWeakInteraction>>(m_sub, "NaivWeakInteraction");


    py::class_<SecondariesCalculator, std::shared_ptr<SecondariesCalculator>>(m_sub, "SecondariesCalculator")
        .def(py::init<std::vector<InteractionType>, ParticleDef const&, Medium const&>())
        .def("n_rnd", &SecondariesCalculator::RequiredRandomNumbers)
        .def("calculate_secondaries", &SecondariesCalculator::CalculateSecondaries);


    m_sub.def("make_secondary",
            [](InteractionType type, ParticleDef const& p, Medium const& m) {
                return std::shared_ptr<secondaries::Parametrization>(secondaries::DefaultFactory::Create(type, p, m));
            });
}
