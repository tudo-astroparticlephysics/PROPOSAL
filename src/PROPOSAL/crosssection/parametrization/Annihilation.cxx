
#include <cassert>
#include <cmath>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crosssection/parametrization/Annihilation.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/Particle.h"

using std::make_tuple;
using namespace PROPOSAL;

double crosssection::Annihilation::GetLowerEnergyLim(
    ParticleDef const& p_def) const
{
    return p_def.mass * 2.f;
}

KinematicLimits crosssection::Annihilation::GetKinematicLimits(
    ParticleDef const& p_def, Component const&, double energy) const
{
    // Limits according to simple 2->2 body interactions
    assert(energy >= p_def.mass);
    auto gamma = energy / p_def.mass;
    auto aux = std::sqrt((gamma - 1.) / (gamma + 1.));
    auto kin_lim = KineamticLimits();
    kin_lim.v_min = 0.5 * (1. - aux);
    kin_lim.v_max = 0.5 * (1. + aux);
    return kin_lim;
}

double crosssection::AnnihilationHeitler::DifferentialCrossSection(
    ParticleDef const& p_def, Component const& comp, double energy,
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
