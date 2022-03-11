
#include <cmath>

#include "PROPOSAL/crosssection/parametrization/PhotoMuPairProduction.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/particle/Particle.h"

using namespace PROPOSAL;
using std::make_tuple;

double crosssection::PhotoMuPairProduction::GetLowerEnergyLim(
        const ParticleDef&) const noexcept {
    return 2. * MMU * SQRTE;
}

crosssection::KinematicLimits crosssection::PhotoMuPairProduction::GetKinematicLimits(
        const ParticleDef&, const Component&, double energy) const noexcept {
    // x is the integration variable here
    KinematicLimits lim{};

    if (energy <= 2. * MMU * SQRTE) {
        lim.v_min = 0.5;
        lim.v_max = 0.5;
    } else {
        lim.v_min = 0.5 * (1. - std::sqrt(1. - 2 * SQRTE * MMU / energy));
        lim.v_max = 0.5 * (1. + std::sqrt(1. - 2 * SQRTE * MMU / energy));
    }
    return lim;
}

std::unique_ptr<crosssection::Parametrization<Component>>
crosssection::PhotoMuPairBurkhardtKelnerKokoulin::clone() const
{
    using param_t = std::remove_cv_t<std::remove_pointer_t<decltype(this)>>;
    return std::make_unique<param_t>(*this);
}

// Parmatrization similar to CERN–SL–2002–016 (https://cds.cern.ch/record/558831/files/sl-2002-016.pdf)
// with added atomic form factors from (Phys. Atom. Nucl. 60 (1997), 576)
double crosssection::PhotoMuPairBurkhardtKelnerKokoulin::DifferentialCrossSection(
        const ParticleDef&, const Component& comp, double energy, double x) const
{
    double aux1, aux2;

    auto delta = MMU * MMU / (2. * energy * x * (1. - x));
    auto B = comp.GetLogConstant();
    auto Bprime = comp.GetBPrime();
    auto Z = comp.GetNucCharge();
    auto Z3 = std::pow(Z, -1. / 3);
    auto D = 1.54 * std::pow(comp.GetAtomicNum(), 0.27);

    // calculate elastic formfactor Phi
    double Phi;
    aux1 = MMU / ME * B * Z3 / (1. + B * Z3 * SQRTE * delta / ME);
    aux2 = D / (1. + delta * (D * SQRTE - 2.) / MMU);
    // influence of the nucleus excitation
    if (Z > 1)
        Phi = std::log(aux1) - (1. - 1. / Z) * std::log(aux2);
    else
        Phi = std::log(aux1) - std::log(aux2);

    // calculate contribution of atomic electrons
    aux1 = (MMU / delta) / (delta * MMU / (ME * ME) + SQRTE);
    aux2 = 1. + 1. / (delta * SQRTE * Bprime * Z3 * Z3 / ME);
    auto electron_contribution = 1. / Z * (std::log(aux1) - std::log(aux2));

    auto aux = x * x + std::pow(1. - x, 2.) + 2./3. * x * (1. - x);
    auto result = 4. * Z * Z * ALPHA * std::pow(RE * ME / MMU, 2.)
        * (Phi + electron_contribution) * aux;

    return std::max(NA / comp.GetAtomicNum() * result, 0.);
}
