
#include <cassert>
#include <cmath>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/parametrization/Annihilation.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/Particle.h"

using std::make_tuple;
using namespace PROPOSAL;

crosssection::Annihilation::Annihilation()
    : Parametrization(InteractionType::Annihilation, "annihililation")
{
}

double crosssection::Annihilation::GetLowerEnergyLim(
    const ParticleDef& p_def) const noexcept
{
    return p_def.mass * 2.f;
}

tuple<double, double> crosssection::Annihilation::GetKinematicLimits(
    const ParticleDef& p_def, const Component&, double energy) const noexcept
{
    // Limits according to simple 2->2 body interactions
    assert(energy >= p_def.mass);
    auto gamma = energy / p_def.mass;
    auto aux = std::sqrt((gamma - 1.) / (gamma + 1.));
    auto vmin = 0.5 * (1. - aux);
    auto vmax = 0.5 * (1. + aux);
    return make_tuple(vmin, vmax);
}

double crosssection::AnnihilationHeitler::DifferentialCrossSection(
    const ParticleDef& p_def, const Component& comp, double energy,
    double v) const
{
    // W. Heitler. The Quantum Theory of Radiation, Clarendon Press, Oxford
    // (1954) Adapted from Geant4 PhysicsReferenceManual

    // v = energy of photon1 / total available energy
    // with the total available energy being the sum of the total positron
    // energy and the electron mass

    assert(energy >= p_def.mass);
    auto gamma = energy / p_def.mass;
    auto aux = 1. + (2. * gamma) / std::pow(gamma + 1., 2.) - v
        - 1. / std::pow(gamma + 1., 2.) * 1. / v;
    aux *= NA * comp.GetNucCharge() / comp.GetAtomicNum() * PI * RE * RE
        / (gamma - 1.) * 1. / v; // TODO: prefactors

    return aux;
}
