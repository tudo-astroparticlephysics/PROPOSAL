#include "PROPOSAL/secondaries/weakinteraction/NaivWeakInteraction.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/Particle.h"

#include <cmath>
#include <stdexcept>

using namespace PROPOSAL;

vector<ParticleState>
secondaries::NaivWeakInteraction::CalculateSecondaries(StochasticLoss loss,
                                                       const Component&,
                                                       vector<double>&)
{
    // TODO: Treatment of hadronic parts of interaction
    auto sec = std::vector<ParticleState>();
    sec.emplace_back(static_cast<ParticleType>(weak_partner_type),
                     loss.position, loss.direction, loss.loss_energy, loss.time,
                     0.);
    return sec;
}
