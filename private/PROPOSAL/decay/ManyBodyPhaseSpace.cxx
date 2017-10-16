
// #include "PROPOSAL/Constants.h"

#include <algorithm>    // std::sort

#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/decay/ManyBodyPhaseSpace.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Output.h"

using namespace PROPOSAL;


const std::string ManyBodyPhaseSpace::name_ = "ManyBodyPhaseSpace";

ManyBodyPhaseSpace::ManyBodyPhaseSpace(std::vector<ParticleDef> daughters)
    : DecayChannel()
    , daughters_(daughters)
    , daughter_masses_()
    , sum_daughter_masses_(0.0)
{
    for(std::vector<ParticleDef>::iterator it = daughters_.begin(); it != daughters_.end(); ++it)
    {
        daughter_masses_.push_back(it->mass);
        sum_daughter_masses_ += it->mass;
    }
}

ManyBodyPhaseSpace::~ManyBodyPhaseSpace()
{
}

ManyBodyPhaseSpace::ManyBodyPhaseSpace(const ManyBodyPhaseSpace& mode)
    : DecayChannel(mode)
    , daughters_(mode.daughters_)
    , daughter_masses_(mode.daughter_masses_)
    , sum_daughter_masses_(mode.sum_daughter_masses_)
{
}

bool ManyBodyPhaseSpace::compare(const DecayChannel& channel) const
{
    const ManyBodyPhaseSpace* many_body = dynamic_cast<const ManyBodyPhaseSpace*>(&channel);

    if (!many_body)
    {
        return false;
    }
    else if (daughter_masses_.size() != many_body->daughter_masses_.size())
    {
        return false;
    }
    else if (daughter_masses_ != many_body->daughter_masses_)
    {
        return false;
    }
    else
    {
        return true;
    }
}

DecayChannel::DecayProducts ManyBodyPhaseSpace::Decay(Particle& particle)
{
    double parent_mass = particle.GetMass();
    std::vector<double> virtual_masses;

    std::vector<double> randoms;
    randoms.reserve(daughters_.size());

    DecayProducts products;
    products.reserve(daughters_.size());

    for (std::vector<ParticleDef>::const_iterator iter = daughters_.begin(); iter != daughters_.end(); ++iter)
    {
        products.push_back(new Particle(*iter));
    }

    // Create sorted random numbers
    randoms.push_back(0.0);

    for (unsigned int i = 0; i < daughter_masses_.size() - 2; ++i)
    {
        randoms.push_back(RandomGenerator::Get().RandomDouble());
    }

    randoms.push_back(1.0);

    std::sort(randoms.begin(), randoms.end());

    // Calculate virtual masses
    double intermediate_mass = 0.0;
    for (unsigned int i = 0; i < daughter_masses_.size(); ++i)
    {
        intermediate_mass += daughter_masses_[i];
        virtual_masses.push_back(intermediate_mass + randoms[i] * (parent_mass - sum_daughter_masses_));
    }

    // Calculate first momentum R2
    unsigned int i = 1;

    double momentum  = Momentum(virtual_masses[i], virtual_masses[i - 1], daughter_masses_[i]);
    Vector3D direction = GenerateRandomDirection();

    products[i]->SetDirection(direction);
    products[i]->SetMomentum(momentum);

    products[i - 1]->SetDirection(-direction);
    products[i - 1]->SetMomentum(momentum);

    for (unsigned int i = 2; i < daughter_masses_.size(); ++i)
    {
        // Generate momentum

        double momentum  = Momentum(virtual_masses[i], virtual_masses[i - 1], daughter_masses_[i]);

        products[i]->SetDirection(GenerateRandomDirection());
        products[i]->SetMomentum(momentum);

        // Boost previous particles to new frame

        double M    = virtual_masses[i - 1];
        double beta = momentum / std::sqrt(momentum * momentum + M * M);

        for (unsigned int s = 0; s < i; ++s)
        {
            Boost(*products[s], products[i]->GetDirection(), beta);
        }
    }

    // Boost all daughters to parent frame
    Boost(products, particle.GetDirection(), particle.GetMomentum() / particle.GetEnergy());

    CopyParticleProperties(products, particle);

    return products;
}

// ------------------------------------------------------------------------- //
// Print
// ------------------------------------------------------------------------- //

void ManyBodyPhaseSpace::print(std::ostream& os) const
{
    os << "Used particle definitions:" << '\n';
    for (std::vector<ParticleDef>::const_iterator it = daughters_.begin(); it != daughters_.end(); ++it)
    {
        os << it->name << '\t' << it->mass << '\n';
    }
    os << "Sum of masses: " << sum_daughter_masses_ << '\n';
}

/******************************************************************************
*                                  Builder                                    *
******************************************************************************/

ManyBodyPhaseSpace::Builder::Builder()
    : daughters_()
{
}

ManyBodyPhaseSpace::Builder::Builder(const Builder& builder)
    :daughters_(builder.daughters_)
{
}

ManyBodyPhaseSpace::Builder::~Builder()
{
}

// ------------------------------------------------------------------------- //
ManyBodyPhaseSpace::Builder& ManyBodyPhaseSpace::Builder::addDaughter(const ParticleDef& daughter)
{
    daughters_.push_back(daughter);
    return *this;
}

// ------------------------------------------------------------------------- //
ManyBodyPhaseSpace ManyBodyPhaseSpace::Builder::build()
{
    return ManyBodyPhaseSpace(daughters_);
}
