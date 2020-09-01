#include "PROPOSAL/secondaries/weakinteraction/NaivWeakInteraction.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/Particle.h"

#include <cmath>
#include <stdexcept>

using namespace PROPOSAL;

vector<Loss::secondary_t>
secondaries::NaivWeakInteraction::CalculateSecondaries(double,
    Loss::secondary_t loss, const Component&,
    vector<double>)
{
    std::get<Loss::TYPE>(loss) = weak_partner_type;
    auto sec = vector<Loss::secondary_t>{ move(loss) };
    return sec;
}
