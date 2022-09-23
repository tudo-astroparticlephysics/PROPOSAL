
#include <cmath>
#include <memory>
#include <type_traits>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crosssection/parametrization/Annihilation.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/medium/Medium.h"

using std::make_tuple;
using namespace PROPOSAL;

InteractionType crosssection::Annihilation::GetInteractionType() const noexcept {
    return InteractionType::Annihilation;
}

size_t crosssection::Annihilation::GetHash(const ParticleDef& p_def, const Medium &m, cut_ptr) const noexcept {
    auto combined_hash = m.GetHash();
    hash_combine(combined_hash, p_def.mass, hash);
    return combined_hash;
}

double crosssection::Annihilation::GetLowerEnergyLim(const ParticleDef& p_def, const Medium&, cut_ptr) const {
    return p_def.mass;
};

crosssection::AnnihilationHeitler::AnnihilationHeitler() {
    hash_combine(hash, std::string(crosssection::ParametrizationName<AnnihilationHeitler>::value));
};

double crosssection::AnnihilationHeitler::CalculatedNdx(double energy, size_t comp_hash, const ParticleDef& p_def,
                                                              const Medium& medium, cut_ptr cut) {
    // integrated form of Heitler Annihilation, cf. Geant4 PhysicsReferenceManual
    if (energy <= Annihilation::GetLowerEnergyLim(p_def, medium, cut))
        return 0.;

    auto comp = Component::GetComponentForHash(comp_hash);
    auto gamma = energy / p_def.mass;
    auto weight = detail::weight_component(medium, comp);
    auto aux = (gamma * gamma + 4 * gamma + 1) / (gamma * gamma - 1)
            * std::log(gamma + std::sqrt(gamma * gamma - 1)) - (gamma + 3) / std::sqrt(gamma * gamma - 1);
    aux *= NA * comp.GetNucCharge() / comp.GetAtomicNum() * PI * RE * RE / (gamma + 1.) / weight;

    return aux;
}

double crosssection::AnnihilationHeitler::CalculatedNdx(double energy, const ParticleDef& p, const Medium& m, cut_ptr cut) {
    double sum = 0.;
    for (auto& rate : CalculatedNdx_PerTarget(energy, p, m, cut))
        sum += rate.second;
    return sum;
}

std::vector<std::pair<size_t, double>> crosssection::AnnihilationHeitler::CalculatedNdx_PerTarget(
        double energy, const ParticleDef& p, const Medium& m, cut_ptr cut) {
    std::vector<std::pair<size_t, double>> rates = {};
    for (auto& comp : m.GetComponents()) {
        double rate = CalculatedNdx(energy, comp.GetHash(), p, m, cut);
        rates.push_back({comp.GetHash(), rate});
    }
    return rates;
}

std::unique_ptr<crosssection::ParametrizationDirect> crosssection::AnnihilationHeitler::clone() const {
    using param_t = std::remove_cv_t<std::remove_pointer_t<decltype(this)>>;
    return std::make_unique<param_t>(*this);
}
