
#include <cassert>
#include <cmath>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/parametrization/Annihilation.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/Particle.h"

using namespace PROPOSAL;

Annihilation::Annihilation()
    : Parametrization(InteractionType::Annihilation, "annihililation")
{
}

Parametrization::KinematicLimits Annihilation::GetKinematicLimits(
    const ParticleDef& p_def, const Component& comp, double energy)
{
    // Limits according to simple 2->2 body interactions

    assert(energy >= p_def.mass);
    auto gamma = energy / p_def.mass;
    auto aux = std::sqrt((gamma - 1.) / (gamma + 1.));

    auto vmin = 0.5 * (1. - aux);
    auto vmax = 0.5 * (1. + aux);

    return KinematicLimits(vmin, vmax);
}

double AnnihilationHeitler::DifferentialCrossSection(
    const ParticleDef& p_def, const Component& comp, double energy, double v)
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
