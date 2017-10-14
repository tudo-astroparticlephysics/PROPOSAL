
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/decay/TwoBodyPhaseSpace.h"
#include "PROPOSAL/particle/Particle.h"

using namespace PROPOSAL;

const std::string TwoBodyPhaseSpace::name_ = "TwoBodyPhaseSpace";

TwoBodyPhaseSpace::TwoBodyPhaseSpace(double m1, double m2)
    : DecayChannel()
    , first_daughter_mass_(m1)
    , second_daughter_mass_(m2)
{
}

TwoBodyPhaseSpace::~TwoBodyPhaseSpace()
{
}

TwoBodyPhaseSpace::TwoBodyPhaseSpace(const TwoBodyPhaseSpace& mode)
    : DecayChannel(mode)
    , first_daughter_mass_(mode.first_daughter_mass_)
    , second_daughter_mass_(mode.second_daughter_mass_)
{
}

bool TwoBodyPhaseSpace::compare(const DecayChannel& channel) const
{
    const TwoBodyPhaseSpace* two_body = dynamic_cast<const TwoBodyPhaseSpace*>(&channel);

    if (!two_body)
        return false;
    else if (first_daughter_mass_ != two_body->first_daughter_mass_)
        return false;
    else if (second_daughter_mass_ != two_body->second_daughter_mass_)
        return false;
    else
        return true;
}

DecayChannel::DecayProducts TwoBodyPhaseSpace::Decay(Particle* particle)
{
    double parent_mass = particle->GetMass();
    double el = (parent_mass * parent_mass + first_daughter_mass_ * first_daughter_mass_) / (2 * parent_mass);

    double arnd = RandomGenerator::Get().RandomDouble();

    double final_energy = el * (particle->GetEnergy() / particle->GetMass()) +
         sqrt(el * el - first_daughter_mass_ * first_daughter_mass_) * (particle->GetMomentum() / particle->GetMass()) * (2 * arnd - 1);

    // return el;
    // Create products
    Particle* product_particle = new Particle(ParticleDef::Builder().SetEMinus().build());
    product_particle->SetEnergy(final_energy);
    product_particle->SetPosition(particle->GetPosition());
    product_particle->SetDirection(particle->GetDirection());
    product_particle->SetParticleId(particle->GetParticleId() + 1);
    product_particle->SetParentParticleId(particle->GetParentParticleId());
    product_particle->SetTime(particle->GetTime());
    product_particle->SetParentParticleEnergy(particle->GetEnergy());

    DecayProducts decay_products;
    decay_products.push_back(product_particle);

    return decay_products;
}

// ------------------------------------------------------------------------- //
// Print
// ------------------------------------------------------------------------- //

void TwoBodyPhaseSpace::print(std::ostream& os) const
{
    os << "Used masses:" << '\n';
    os << "First daughter mass: " << first_daughter_mass_ << '\n';
    os << "Second daughter mass: " << second_daughter_mass_ << '\n';
}
