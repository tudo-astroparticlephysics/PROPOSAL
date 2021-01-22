#include "PROPOSAL/secondaries/parametrization/weakinteraction/NaivWeakInteraction.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/Particle.h"

#include <cmath>
#include <stdexcept>

using namespace PROPOSAL;

std::vector<ParticleState>
secondaries::NaivWeakInteraction::CalculateSecondaries(StochasticLoss loss,
                                                       const Component&,
                                                       std::vector<double>&)
{
    // TODO: Treatment of hadronic parts of interaction
    auto sec = std::vector<ParticleState>();
    sec.emplace_back(static_cast<ParticleType>(weak_partner_type),
                     loss.position, loss.direction, loss.energy, loss.time,
                     0.);
    return sec;
}
