
#include "PROPOSAL/decay/TwoBodyPhaseSpace.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/particle/Particle.h"

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

DecayChannel::DecayProducts TwoBodyPhaseSpace::Decay(const Particle& particle)
{
    DecayProducts products;
    products.push_back(new Particle(first_daughter_));
    products.push_back(new Particle(second_daughter_));

    double momentum    = Momentum(particle.GetMass(), first_daughter_.mass, second_daughter_.mass);
    Vector3D direction = GenerateRandomDirection();

    products[0]->SetDirection(direction);
    products[0]->SetMomentum(momentum);

    products[1]->SetDirection(-direction);
    products[1]->SetMomentum(momentum);

    Boost(products, particle.GetDirection(), particle.GetMomentum() / particle.GetEnergy());

    CopyParticleProperties(products, particle);

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
