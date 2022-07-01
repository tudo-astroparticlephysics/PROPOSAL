#include "PROPOSAL/secondaries/parametrization/photoeffect/NaivPhotoeffect.h"
#include "PROPOSAL/Constants.h"

#include <cmath>

using namespace PROPOSAL;

std::vector<ParticleState> secondaries::NaivPhotoeffect::CalculateSecondaries(
        StochasticLoss loss, const Component& comp, std::vector<double>&)
{
    double I = pow(comp.GetNucCharge() * ALPHA, 2) * ME / 2; // subtract ionization energy of K-shell electron
    double E = ME + loss.parent_particle_energy - I;

    auto sec = std::vector<ParticleState>();
    // TODO: Do we want to do something better here?
    // TODO: Do we need to make sure that the created electron has a valid energy that is above its rest mass? Or is this already secured?
    // assume electron inherits direction of photon
    sec.emplace_back(ParticleType::EMinus, loss.position, loss.direction, E, loss.time, 0.);
    return sec;
}
