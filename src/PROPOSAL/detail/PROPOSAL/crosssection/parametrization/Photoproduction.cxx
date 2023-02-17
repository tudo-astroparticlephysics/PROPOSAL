#include "PROPOSAL/crosssection/parametrization/Photoproduction.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/math/Interpolant.h"

using namespace PROPOSAL;

InteractionType crosssection::Photoproduction::GetInteractionType() const noexcept {
    return InteractionType::Photoproduction;
}

size_t crosssection::Photoproduction::GetHash(const ParticleDef&, const Medium &m, cut_ptr) const noexcept {
    auto combined_hash = m.GetHash();
    hash_combine(combined_hash, hash);
    return combined_hash;
}

double crosssection::Photoproduction::GetLowerEnergyLim(const ParticleDef&, const Medium& medium, cut_ptr) const {
    double min_cut_off = INF;
    for (const auto& c : medium.GetComponents())
        min_cut_off = std::min(min_cut_off, GetCutOff(c));
    return min_cut_off;
};

double crosssection::Photoproduction::GetCutOff(const Component& comp) const {
    // pion production threshold
    auto m_N = comp.GetAverageNucleonWeight();
    return MPI + MPI * MPI / (2. * m_N);
}

double crosssection::Photoproduction::PhotonAtomCrossSection(double energy, const Component& comp) {
    auto cross_photon_nucleon = PhotonNucleonCrossSection(energy, comp);
    if (cross_photon_nucleon == 0.)
        return 0.;
    return comp.GetAtomicNum() * cross_photon_nucleon * (0.75 * ShadowingFactor(energy, comp) + 0.25);
}

double crosssection::Photoproduction::ShadowingFactor(double energy, const Component& comp) {
    if (comp.GetNucCharge() == 1) {
        return 1.;
    } else {
        double x = 0.00282 * std::pow(comp.GetAtomicNum(), 1. / 3) * PhotonNucleonCrossSection(energy, comp);
        return 3 / (x * x * x) * ( (x * x) / 2 - 1 + std::exp(-x) * (1 + x));
    }
}

double crosssection::Photoproduction::CalculatedNdx(
        double energy, size_t comp_hash, const ParticleDef&, const Medium& m, cut_ptr) {
        auto comp = Component::GetComponentForHash(comp_hash);
        auto weight = detail::weight_component(m, comp);
        return NA / comp.GetAtomicNum() * 1e-30 * PhotonAtomCrossSection(energy, comp) / weight;
}

double crosssection::Photoproduction::CalculatedNdx(double energy, const ParticleDef& p, const Medium& m, cut_ptr cut) {
    double sum = 0.;
    for (auto& rate : CalculatedNdx_PerTarget(energy, p, m, cut))
        sum += rate.second;
    return sum;
}

std::vector<std::pair<size_t, double>> crosssection::Photoproduction::CalculatedNdx_PerTarget(
        double energy, const ParticleDef& p, const Medium& m, cut_ptr) {
    std::vector<std::pair<size_t, double>> rates = {};

    for (auto& comp : m.GetComponents()) {
        auto weight = detail::weight_component(m, comp);
        double rate = NA / comp.GetAtomicNum() * 1e-30 * PhotonAtomCrossSection(energy, comp) / weight;
        rates.push_back({comp.GetHash(), rate});
    }
    return rates;
}

// Zeus
crosssection::PhotoproductionZeus::PhotoproductionZeus() {
    hash_combine(hash, std::string(crosssection::ParametrizationName<PhotoproductionZeus>::value));
};

std::unique_ptr<crosssection::ParametrizationDirect> crosssection::PhotoproductionZeus::clone() const {
    using param_t = std::remove_cv_t<std::remove_pointer_t<decltype(this)>>;
    return std::make_unique<param_t>(*this);
}

// ------------------------------------------------------------------------- //
// Zeus Collaboration, Breitweg et al
// Eur. Phys. J. C 7 (1999), 609
// eq. 6
// ------------------------------------------------------------------------- //
double crosssection::PhotoproductionZeus::PhotonNucleonCrossSection(double energy, const Component& comp) {
    if (energy < GetCutOff(comp))
        return 0.;
    double aux;
    auto nu = energy * 1e-3; // from MeV to GeV
    aux = 2e-3 * nu * comp.GetAverageNucleonWeight(); // NucleonWeight from MeV to GeV
    aux = 63.5 * std::pow(aux, 0.097) + 145 * std::pow(aux, -0.5);
    return aux; // return value in μb
}

// BezrukovBugaev
crosssection::PhotoproductionBezrukovBugaev::PhotoproductionBezrukovBugaev() {
    hash_combine(hash, std::string(crosssection::ParametrizationName<PhotoproductionBezrukovBugaev>::value));
};

std::unique_ptr<crosssection::ParametrizationDirect> crosssection::PhotoproductionBezrukovBugaev::clone() const {
    using param_t = std::remove_cv_t<std::remove_pointer_t<decltype(this)>>;
    return std::make_unique<param_t>(*this);
}

// ------------------------------------------------------------------------- //
// Bezrukov, Bugaev
// Sov. J. Nucl. Phys. 32 (1980), 847
// eq. 21
// ------------------------------------------------------------------------- //
double crosssection::PhotoproductionBezrukovBugaev::PhotonNucleonCrossSection(double energy, const Component& comp) {
    if (energy < GetCutOff(comp))
        return 0.;
    auto nu = energy * 1e-3; // from MeV to GeV
    double aux;
    aux = std::log(0.0213 * nu);
    aux = 114.3 + 1.647 * aux * aux;
    return aux; // return value in μb
}

// Caldwell
crosssection::PhotoproductionCaldwell::PhotoproductionCaldwell() {
    hash_combine(hash, std::string(crosssection::ParametrizationName<PhotoproductionCaldwell>::value));
};

std::unique_ptr<crosssection::ParametrizationDirect> crosssection::PhotoproductionCaldwell::clone() const {
    using param_t = std::remove_cv_t<std::remove_pointer_t<decltype(this)>>;
    return std::make_unique<param_t>(*this);
}

// ------------------------------------------------------------------------- //
// fit of Caldwell at al.
// Phys. Rev Let. 42 (1979), 553
// Table 1
// ------------------------------------------------------------------------- //
double crosssection::PhotoproductionCaldwell::PhotonNucleonCrossSection(double energy, const Component& comp) {
    if (energy < GetCutOff(comp))
        return 0.;
    auto nu = energy * 1e-3; // from MeV to GeV
    return 49.2 + 11.1 * std::log(nu) + 151.8 / std::sqrt(nu); // return value in μb
}

// Kokoulin
crosssection::PhotoproductionKokoulin::PhotoproductionKokoulin() {
    hash_combine(hash, std::string(crosssection::ParametrizationName<PhotoproductionKokoulin>::value));
};

std::unique_ptr<crosssection::ParametrizationDirect> crosssection::PhotoproductionKokoulin::clone() const {
    using param_t = std::remove_cv_t<std::remove_pointer_t<decltype(this)>>;
    return std::make_unique<param_t>(*this);
}

double crosssection::PhotoproductionKokoulin::PhotonNucleonCrossSection(double energy, const Component& comp) {
    if (energy < GetCutOff(comp))
        return 0.;
    auto nu = energy * 1e-3; // from MeV to GeV
    if (nu <= 200.) {
        if (nu <= 17.) {
            // Bezrukov, Bugaev, Sov. J. Nucl. Phys. 33 (1981), 635
            return 96.1 + 82. / std::sqrt(nu);
        } else {
            return PhotoproductionBezrukovBugaev::PhotonNucleonCrossSection(energy, comp);
        }
    } else {
        return PhotoproductionCaldwell::PhotonNucleonCrossSection(energy, comp);
    }
}

// Rhode

crosssection::PhotoproductionRhode::PhotoproductionRhode() {
        std::vector<double> x = { 0, 0.1, 0.144544, 0.20893, 0.301995, 0.436516,
        0.630957, 0.912011, 1.31826, 1.90546, 2.75423, 3.98107, 5.7544, 8.31764,
        12.0226, 17.378, 25.1189, 36.3078, 52.4807, 75.8577, 109.648, 158.489,
        229.087, 331.131, 478.63, 691.831, 1000, 1445.44, 2089.3, 3019.95,
        4365.16, 6309.58, 9120.12, 13182.6, 19054.6, 27542.3, 39810.8, 57544,
        83176.4, 120226, 173780, 251188, 363078, 524807, 758576, 1.09648e+06,
        1.58489e+06, 2.29086e+06, 3.3113e+06, 4.78628e+06, 6.91828e+06,
        9.99996e+06 };

    std::vector<double> y = { 0, 0.0666667, 0.0963626, 159.74, 508.103, 215.77,
        236.403, 201.919, 151.381, 145.407, 132.096, 128.546, 125.046, 121.863,
        119.16, 117.022, 115.496, 114.607, 114.368, 114.786, 115.864, 117.606,
        120.011, 123.08, 126.815, 131.214, 136.278, 142.007, 148.401, 155.46,
        163.185, 171.574, 180.628, 190.348, 200.732, 211.782, 223.497, 235.876,
        248.921, 262.631, 277.006, 292.046, 307.751, 324.121, 341.157, 358.857,
        377.222, 396.253, 415.948, 436.309, 457.334, 479.025 };

    interpolant_ = std::make_shared<Interpolant>(x, y, 4, false, false);
    hash_combine(hash, std::string(crosssection::ParametrizationName<PhotoproductionRhode>::value));
}

std::unique_ptr<crosssection::ParametrizationDirect> crosssection::PhotoproductionRhode::clone() const {
    using param_t = std::remove_cv_t<std::remove_pointer_t<decltype(this)>>;
    return std::make_unique<param_t>(*this);
}

double crosssection::PhotoproductionRhode::PhotonNucleonCrossSection(double energy, const Component& comp) {
    if (energy < GetCutOff(comp))
        return 0.;
    auto nu = energy * 1e-3; // from MeV to GeV
    if (nu <= 0.1) {
        return 0.; // do not extrapolate
    } else if (nu <= 200.) {
        return std::max(interpolant_->InterpolateArray(nu), 0.);
    } else {
        return PhotoproductionCaldwell::PhotonNucleonCrossSection(energy, comp);
    }
}

// Heck
crosssection::PhotoproductionHeck::PhotoproductionHeck() : resonances_M{{1.231, 1.515, 1.680}},
                                                           resonances_Gamma{{0.11, 0.11, 0.125}},
                                                           resonances_Sigma0{{31.125, 25.567, 17.508}},
                                                           resonances_w{{0.17, 0.38, 0.38}} {
    hash_combine(hash, std::string(crosssection::ParametrizationName<PhotoproductionHeck>::value));
};

std::unique_ptr<crosssection::ParametrizationDirect> crosssection::PhotoproductionHeck::clone() const {
    using param_t = std::remove_cv_t<std::remove_pointer_t<decltype(this)>>;
    return std::make_unique<param_t>(*this);
}

// ------------------------------------------------------------------------- //
// Formulae of continuum are from https://web.iap.kit.edu/heck/publications/kit-swp-6_vectormeson.pdf
// The resonances are implemented according to https://arxiv.org/pdf/astro-ph/9903478.pdf using the
// Delta(1232), N(1520) and N(1680) resonance.
// To describe shadowing, the parametrization by BezrukovBugaev is used
// ------------------------------------------------------------------------- //
double crosssection::PhotoproductionHeck::PhotonNucleonCrossSection(double energy, const Component& comp) {
    if (energy < GetCutOff(comp))
        return 0.;
    energy *= 1e-3; // energy in GeV
    double m_n = comp.GetAverageNucleonWeight() * 1e-3; // in GeV
    double s = m_n * m_n + 2 * m_n * energy; // calculate squared center of mass energy in GeV^2
    double s0 = std::pow(1.0761, 2); // pion production threshold

    double total_cross;

    // calculate continuum of crosssection
    if (std::sqrt(s) < 19.39)
        total_cross = (73.7 * std::pow(s, 0.073) + 191.7 * std::pow(s, -0.602)) * std::sqrt(1 - s0/s); // in µb
    else
        total_cross = 59.3 * std::pow(s, 0.093) + 120.2 * std::pow(s, -0.358); // in µb

    // add resonances
    for (unsigned int i = 0; i<resonances_M.size(); i++) {
        total_cross += sigma_bw(s, energy, resonances_M.at(i), resonances_Gamma.at(i), resonances_Sigma0.at(i))
                * Q(energy, 0.152, resonances_w.at(i));
    }

    return total_cross;
}

double crosssection::PhotoproductionHeck::sigma_bw(double s, double E, double M, double Gamma, double Sigma0) {
        return s / (E * E) * Sigma0 * Gamma * Gamma * s / (std::pow(s - M*M, 2) + Gamma * Gamma * s);
}

double crosssection::PhotoproductionHeck::Q(double E, double E_th, double w) {
    if (E < E_th)
        return 0;
    if (E >= w + E_th)
        return 1;
    return (E - E_th) / w;
}


// HeckC7Shadowing
crosssection::PhotoproductionHeckC7Shadowing::PhotoproductionHeckC7Shadowing() {
    hash_combine(hash, std::string(crosssection::ParametrizationName<PhotoproductionHeckC7Shadowing>::value));
};

std::unique_ptr<crosssection::ParametrizationDirect> crosssection::PhotoproductionHeckC7Shadowing::clone() const {
    using param_t = std::remove_cv_t<std::remove_pointer_t<decltype(this)>>;
    return std::make_unique<param_t>(*this);
}

// ------------------------------------------------------------------------- //
// Use simplified shadowing method according to https://web.iap.kit.edu/heck/publications/kit-swp-6_vectormeson.pdf
// or Phys. Rev. D 32 (1985) 1244 (https://journals.aps.org/prd/pdf/10.1103/PhysRevD.32.1244)
// ------------------------------------------------------------------------- //
double crosssection::PhotoproductionHeckC7Shadowing::PhotonAtomCrossSection(double energy, const Component& comp) {
    auto cross_photon_nucleon = PhotonNucleonCrossSection(energy, comp);
    if (cross_photon_nucleon == 0.)
        return 0.;
    return cross_photon_nucleon * std::pow(comp.GetAtomicNum(), 0.91);
};