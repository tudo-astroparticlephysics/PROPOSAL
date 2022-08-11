#include "PROPOSAL/secondaries/parametrization/photoeffect/PhotoeffectNoDeflection.h"
#include "PROPOSAL/Constants.h"

#include <cmath>
#include <cassert>

using namespace PROPOSAL;

std::vector<ParticleState> secondaries::PhotoeffectNoDeflection::CalculateSecondaries(
        StochasticLoss loss, const Component& comp, std::vector<double>&)
{
    double I = pow(comp.GetNucCharge() * ALPHA, 2) * ME / 2; // subtract ionization energy of K-shell electron
    double E = ME + loss.parent_particle_energy - I;
    assert(E >= ME);

    auto sec = std::vector<ParticleState>();
    // assume electron inherits direction of photon (this is also done in standard EGS4)
    sec.emplace_back(ParticleType::EMinus, loss.position, loss.direction, E, loss.time, 0.);
    return sec;
}
