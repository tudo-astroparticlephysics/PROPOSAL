
#include <cassert>
#include <cmath>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/parametrization/Annihilation.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/particle/Particle.h"

using namespace PROPOSAL;

Annihilation::Annihilation(const ParticleDef& p_def, const component_list& comp)
    : Parametrization(InteractionType::Annihilation, "annihililation", p_def, comp, 2 * ME)
    , za_(calculate_proton_massnumber_fraction(comp))
{
}

Parametrization::KinematicLimits Annihilation::GetKinematicLimits(double energy)
{
    // Limits according to simple 2->2 body interactions

    assert(energy >= particle_mass_);
    auto gamma = energy / particle_mass_;
    auto aux = std::sqrt((gamma - 1.) / (gamma + 1.));

    auto vmin = 0.5 * (1. - aux);
    auto vmax = 0.5 * (1. + aux);

    return KinematicLimits(vmin, vmax);
}

AnnihilationHeitler::AnnihilationHeitler(
    const ParticleDef& p_def, const component_list& comp)
    : Annihilation(p_def, comp)
{
}

double AnnihilationHeitler::DifferentialCrossSection(double energy, double v)
{
    // W. Heitler. The Quantum Theory of Radiation, Clarendon Press, Oxford
    // (1954) Adapted from Geant4 PhysicsReferenceManual

    // v = energy of photon1 / total available energy
    // with the total available energy being the sum of the total positron
    // energy and the electron mass

    assert(energy >= particle_mass_);
    auto gamma = energy / particle_mass_;
    auto aux = 1. + (2. * gamma) / std::pow(gamma + 1., 2.) - v
        - 1. / std::pow(gamma + 1., 2.) * 1. / v;
    aux *= NA * za_ * PI * RE * RE / (gamma - 1.) * 1. / v; // TODO: prefactors

    return aux;
}
