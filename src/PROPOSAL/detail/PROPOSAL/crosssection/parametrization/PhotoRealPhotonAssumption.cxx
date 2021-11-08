#include <cmath>

#include "PROPOSAL/crosssection/parametrization/PhotoRealPhotonAssumption.h"
#include "PROPOSAL/crosssection/parametrization/Photoproduction.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/particle/ParticleDef.h"

#define PHOTO_PARAM_REAL_IMPL(param)                                           \
    crosssection::Photo##param::Photo##param(bool hard_component)              \
        : crosssection::PhotoRealPhotonAssumption(hard_component,              \
          std::make_shared<crosssection::Photoproduction##param>())            \
    {                                                                          \
        hash_combine(hash, std::string(#param));                               \
    }                                                                          \
                                                                               \
    std::unique_ptr<crosssection::Parametrization<Component>>                  \
        crosssection::Photo##param::clone() const                              \
    {                                                                          \
        using param_t                                                          \
            = std::remove_cv_t<std::remove_pointer_t<decltype(this)>>;         \
        return std::make_unique<param_t>(*this);                               \
    }

using namespace PROPOSAL;

crosssection::PhotoRealPhotonAssumption::PhotoRealPhotonAssumption(
        bool hard_component,
        std::shared_ptr<crosssection::Photoproduction> photon_param)
    : crosssection::Photonuclear()
    , hard_component_(hard_component)
    , photon_param_(std::move(photon_param))
{
    hash_combine(hash, hard_component);
}

// ------------------------------------------------------------------------- //
// Bezrukov, Bugaev
// Sov. J. Nucl. Phys. 33 (1981), 635
// eq. 23
// ------------------------------------------------------------------------- //
double crosssection::PhotoRealPhotonAssumption::DifferentialCrossSection(
    const ParticleDef& p_def, const Component& comp, double energy,
    double v) const
{
    double aux, aum, G, t;

    double sgn = photon_param_->PhotonNucleonCrossSection(v * energy, comp);

    const double m1 = 0.54; // in GeV^2
    const double m2 = 1.80; // in GeV^2

    double kappa = 1 - 2 / v + 2 / (v * v);

    // calculate the shadowing factor
    if (comp.GetNucCharge() == 1) {
        G = 1;
    } else {
        // eq. 18
        double tmp = 0.00282 * std::pow(comp.GetAtomicNum(), 1. / 3) * sgn;
        // eq. 3
        G = (3 / tmp) * (0.5 + ((1 + tmp) * std::exp(-tmp) - 1) / (tmp * tmp));
    }

    // enhanced formula by Bugaev Shelpin
    // Phys. Rev. D 67 (2003), 034027
    // eq. 4.6
    G *= 3;
    aux = v * p_def.mass * 1.e-3;
    t = aux * aux / (1 - v);
    aum = p_def.mass * 1.e-3;
    aum *= aum;
    aux = 2 * aum / t;
    aux = G
            * ((kappa + 4 * aum / m1) * std::log(1 + m1 / t)
                - (kappa * m1) / (m1 + t) - aux)
        + ((kappa + 2 * aum / m2) * std::log(1 + m2 / t) - aux)
        + aux * (G * (m1 - 4 * t) / (m1 + t) + (m2 / t) * std::log(1 + t / m2));

    aux *= ALPHA / (8 * PI) * comp.GetAtomicNum() * v * sgn * 1.e-30;

    // hard component by Bugaev, Montaruli, Shelpin, Sokalski
    // Astrop. Phys. 21 (2004), 491
    // in the appendix
    if (hard_component_)
        aux += comp.GetAtomicNum() * 1.e-30
            * HardComponent(p_def).CalculateHardComponent(energy, v);
    else
        aux += comp.GetAtomicNum() * 1.e-30
            * SoftComponent().CalculateHardComponent(energy, v);

    return NA / comp.GetAtomicNum() * p_def.charge * p_def.charge * aux;
}

double crosssection::PhotoRealPhotonAssumption::CalculateParametrization(
        const Component& comp, double energy) const {
    return photon_param_->PhotonNucleonCrossSection(energy, comp);
}

PHOTO_PARAM_REAL_IMPL(Zeus)
PHOTO_PARAM_REAL_IMPL(BezrukovBugaev)
PHOTO_PARAM_REAL_IMPL(Kokoulin)
PHOTO_PARAM_REAL_IMPL(Rhode)

#undef Q2_PHOTO_PARAM_INTEGRAL_IMPL
