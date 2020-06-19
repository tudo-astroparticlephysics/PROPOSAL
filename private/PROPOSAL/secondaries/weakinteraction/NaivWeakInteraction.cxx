#include "PROPOSAL/secondaries/weakinteraction/NaivWeakInteraction.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/Particle.h"

#include <cmath>
#include <stdexcept>

using namespace PROPOSAL;
secondaries::NaivWeakInteraction::NaivWeakInteraction(
    const ParticleDef& p_def, const Medium&)
    : weak_partner_type(p_def.weak_partner)
{
}

vector<Loss::secondary_t>
secondaries::NaivWeakInteraction::CalculateSecondaries(double,
    Loss::secondary_t loss, const Component&,
    vector<double> rnd)
{
    std::get<Loss::TYPE>(loss) = weak_partner_type;
    auto sec = vector<Loss::secondary_t>{ move(loss) };
    return sec;
}
