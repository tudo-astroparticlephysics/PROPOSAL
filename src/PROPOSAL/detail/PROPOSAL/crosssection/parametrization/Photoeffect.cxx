#include "PROPOSAL/crosssection/parametrization/Photoeffect.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/crosssection/CrossSection.h"

using namespace PROPOSAL;

InteractionType crosssection::Photoeffect::GetInteractionType() const noexcept {
    return InteractionType::Photoeffect;
}

size_t crosssection::Photoeffect::GetHash(const ParticleDef&, const Medium &m, cut_ptr) const noexcept {
    auto combined_hash = m.GetHash();
    hash_combine(combined_hash, hash);
    return combined_hash;
}

double crosssection::Photoeffect::GetLowerEnergyLim(const ParticleDef&, const Medium& medium, cut_ptr) const {
    double min_cut_off = INF;
    for (const auto& c : medium.GetComponents())
        min_cut_off = std::min(min_cut_off, GetCutOff(c));
    return min_cut_off;
};

double crosssection::Photoeffect::GetCutOff(const Component& comp) const {
    double I = pow(comp.GetNucCharge() * ALPHA, 2) * ME / 2; // ionization energy of K-shell electron
    return I;
};

double Kshell_Total_Ratio(const Component& comp) {
    // J. H. Hubbell, Tech. Rep. NSRDS-NBS 29, Nat. Bur. Stand. (1969)
    // eq. (2.-2)
    double ln_Z = log(comp.GetNucCharge());
    double res = 1 + 0.01481 * ln_Z * ln_Z - 0.000788 * ln_Z * ln_Z * ln_Z;
    return res;
}

double crosssection::Photoeffect::PhotonAtomCrossSection(double energy, const Component& comp) {
    auto cross_photon_atom = PhotoeffectKshellCrossSection(energy, comp);
    if (cross_photon_atom == 0.)
        return 0.;
    return cross_photon_atom * Kshell_Total_Ratio(comp);
}

double crosssection::Photoeffect::CalculatedNdx(
        double energy, size_t comp_hash, const ParticleDef&, const Medium& m, cut_ptr) {
        auto comp = Component::GetComponentForHash(comp_hash);
        auto weight = detail::weight_component(m, comp);
        return NA / comp.GetAtomicNum() * PhotonAtomCrossSection(energy, comp) / weight;
}

double crosssection::Photoeffect::CalculatedNdx(double energy, const ParticleDef& p, const Medium& m, cut_ptr cut) {
    double sum = 0.;
    for (auto& rate : CalculatedNdx_PerTarget(energy, p, m, cut))
        sum += rate.second;
    return sum;
}

std::vector<std::pair<size_t, double>> crosssection::Photoeffect::CalculatedNdx_PerTarget(
        double energy, const ParticleDef& p, const Medium& m, cut_ptr) {
    std::vector<std::pair<size_t, double>> rates = {};
    for (auto& comp : m.GetComponents()) {
        double rate = 0.;
        if (energy > GetCutOff(comp)) {
            auto weight = detail::weight_component(m, comp);
            rate = NA / comp.GetAtomicNum() * PhotonAtomCrossSection(energy, comp) / weight;
        }
        rates.push_back({comp.GetHash(), rate});
    }
    return rates;
}

// Sauter
crosssection::PhotoeffectSauter::PhotoeffectSauter() {
    hash_combine(hash, std::string(crosssection::ParametrizationName<PhotoeffectSauter>::value));
};

std::unique_ptr<crosssection::ParametrizationDirect> crosssection::PhotoeffectSauter::clone() const {
    using param_t = std::remove_cv_t<std::remove_pointer_t<decltype(this)>>;
    return std::make_unique<param_t>(*this);
}

// ------------------------------------------------------------------------- //
// F. Sauter, Ann. d. Phys. 403 (1931) 454-488
// eq. 40, 41
// cited after W. Heitler, The Quantum Theory of Radiation. 3rd ed.
// Oxford Clarendon Press (1954)
// eq. 17
// ------------------------------------------------------------------------- //
double crosssection::PhotoeffectSauter::PhotoeffectKshellCrossSection(double energy, const Component& comp) {
    const double sigma_T = (8 * PI * RE * RE)/3; // Thomson cross section
    double I = pow(comp.GetNucCharge() * ALPHA, 2) * ME / 2; // ionization energy of K-shell electron
    auto gamma = 1 + (energy - I) / ME;
    auto gamma_beta = sqrt(gamma * gamma - 1);

    auto sigma = 1.5 * sigma_T * pow(comp.GetNucCharge(), 5) * pow(ALPHA, 4) * pow(ME / energy, 5)
        * pow(gamma_beta, 3) * (4/3.0 + gamma * (gamma - 2) / (gamma + 1)
        * (1 - 1 / (gamma * gamma_beta) * log(gamma + gamma_beta)));

    auto Z_alpha_beta = comp.GetNucCharge() * ALPHA * gamma / gamma_beta; // Z*alpha/beta
    auto Z_pi_alpha_beta = PI * Z_alpha_beta;
    // correction factor for non-relativstic electrons
    double non_rel;
    if (Z_pi_alpha_beta < 10.) {
      non_rel = (1 + Z_alpha_beta * Z_alpha_beta) * Z_pi_alpha_beta / sinh(Z_pi_alpha_beta)
        * exp(Z_alpha_beta * (PI - 4 * atan(1 / Z_alpha_beta)));
    } else {
      non_rel = (1 + Z_alpha_beta * Z_alpha_beta) * 2 * Z_pi_alpha_beta
        * exp(Z_alpha_beta * (- 4 * atan(1 / Z_alpha_beta)));
    }
    return sigma * non_rel;
}
