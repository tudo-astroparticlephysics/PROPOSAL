
#include "PROPOSAL/decay/TwoBodyPhaseSpace.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/particle/ParticleDef.h"

using namespace PROPOSAL;

const std::string TwoBodyPhaseSpace::name_ = "TwoBodyPhaseSpace";

TwoBodyPhaseSpace::TwoBodyPhaseSpace(const ParticleDef& first, const ParticleDef& second)
    : DecayChannel()
    , first_daughter_(first)
    , second_daughter_(second)
{
}

TwoBodyPhaseSpace::~TwoBodyPhaseSpace() {}

TwoBodyPhaseSpace::TwoBodyPhaseSpace(const TwoBodyPhaseSpace& mode)
    : DecayChannel(mode)
    , first_daughter_(mode.first_daughter_)
    , second_daughter_(mode.second_daughter_)
{
}

bool TwoBodyPhaseSpace::compare(const DecayChannel& channel) const
{
    const TwoBodyPhaseSpace* two_body = dynamic_cast<const TwoBodyPhaseSpace*>(&channel);

    if (!two_body)
        return false;
    else if (first_daughter_ != two_body->first_daughter_)
        return false;
    else if (second_daughter_ != two_body->second_daughter_)
        return false;
    else
        return true;
}

std::vector<ParticleState> TwoBodyPhaseSpace::Decay(const ParticleDef& p_def, const ParticleState& p_condition)
{
    std::vector<ParticleState> products;
    products.emplace_back((ParticleType)first_daughter_.particle_type, p_condition.position, p_condition.direction, p_condition.energy, p_condition.time, 0);
    products.emplace_back((ParticleType)second_daughter_.particle_type, p_condition.position, p_condition.direction, p_condition.energy, p_condition.time, 0);

    double momentum    = Momentum(p_def.mass, first_daughter_.mass, second_daughter_.mass);
    auto direction = GenerateRandomDirection();

    products[0].direction = direction;
    products[0].SetMomentum(momentum);

    Cartesian3D opposite_direction = -direction;
    products[1].direction = opposite_direction;
    products[1].SetMomentum(momentum);

    // Boost all products in Lab frame (the reason, why the boosting goes in the negative direction of the particle)
    Boost(products, -p_condition.direction, p_condition.energy / p_def.mass, p_condition.GetMomentum() / p_def.mass);

    return products;
}

// ------------------------------------------------------------------------- //
// Print
// ------------------------------------------------------------------------- //

void TwoBodyPhaseSpace::print(std::ostream& os) const
{
    os << "Used masses:" << '\n';
    os << "First daughter:\n" << first_daughter_ << '\n';
    os << "Second daughter:\n" << second_daughter_ << '\n';
}
