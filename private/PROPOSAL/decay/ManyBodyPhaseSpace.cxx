
// #include "PROPOSAL/Constants.h"

#include <algorithm> // std::sort

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/decay/ManyBodyPhaseSpace.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/particle/Particle.h"

#include <boost/math/special_functions/factorials.hpp>

using namespace PROPOSAL;

const std::string ManyBodyPhaseSpace::name_ = "ManyBodyPhaseSpace";

ManyBodyPhaseSpace::ManyBodyPhaseSpace(std::vector<const ParticleDef*> daughters)
    : DecayChannel()
    , daughters_()
    , daughter_masses_()
    , number_of_daughters_(static_cast<int>(daughters.size()))
    , sum_daughter_masses_(0.0)
    , uniform_(true)
{
    for (std::vector<const ParticleDef*>::iterator it = daughters.begin(); it != daughters.end(); ++it)
    {
        daughters_.push_back((*it)->clone());
        daughter_masses_.push_back((*it)->mass);
        sum_daughter_masses_ += (*it)->mass;
    }
}

ManyBodyPhaseSpace::~ManyBodyPhaseSpace()
{
    for (std::vector<const ParticleDef*>::const_iterator it = daughters_.begin(); it != daughters_.end(); ++it)
    {
        delete *it;
    }

    daughters_.clear();
}

ManyBodyPhaseSpace::ManyBodyPhaseSpace(const ManyBodyPhaseSpace& mode)
    : DecayChannel(mode)
    , daughters_()
    , daughter_masses_(mode.daughter_masses_)
    , number_of_daughters_(mode.number_of_daughters_)
    , sum_daughter_masses_(mode.sum_daughter_masses_)
    , uniform_(mode.uniform_)
{
    for (std::vector<const ParticleDef*>::const_iterator it = mode.daughters_.begin(); it != mode.daughters_.end();
         ++it)
    {
        daughters_.push_back((*it)->clone());
    }
}

bool ManyBodyPhaseSpace::compare(const DecayChannel& channel) const
{
    const ManyBodyPhaseSpace* many_body = dynamic_cast<const ManyBodyPhaseSpace*>(&channel);

    if (!many_body)
    {
        return false;
    } else if (daughter_masses_.size() != many_body->daughter_masses_.size())
    {
        return false;
    } else if (number_of_daughters_ != many_body->number_of_daughters_)
    {
        return false;
    } else if (daughter_masses_ != many_body->daughter_masses_)
    {
        return false;
    } else if (uniform_ != many_body->uniform_)
    {
        return false;
    } else
    {
        return true;
    }
}

DecayChannel::DecayProducts ManyBodyPhaseSpace::Decay(const Particle& particle)
{
    double parent_mass = particle.GetMass();

    // Create vector for decay products
    DecayProducts products;
    products.reserve(daughters_.size());

    for (std::vector<const ParticleDef*>::const_iterator it = daughters_.begin(); it != daughters_.end(); ++it)
    {
        products.push_back(new Particle(**it));
    }

    // prefactor for the phase space density
    PhaseSpaceParameters params = GetPhaseSpaceParams(particle.GetParticleDef());

    // precalculated kinematics
    PhaseSpaceKinematics kinematics = CalculateKinematics(params.normalization, parent_mass);

    if (uniform_)
    {
        // do rejection until phase space distribution is uniform
        while (RandomGenerator::Get().RandomDouble() * params.weight_max > kinematics.weight)
        {
            kinematics = CalculateKinematics(params.normalization, parent_mass);
        }
    }

    // Calculate first momentum in R2
    Vector3D direction = GenerateRandomDirection();

    products[1]->SetDirection(direction);
    products[1]->SetMomentum(kinematics.momenta[0]);

    products[0]->SetDirection(-direction);
    products[0]->SetMomentum(kinematics.momenta[0]);

    // Correct the previous momenta
    for (unsigned int i = 2; i < daughter_masses_.size(); ++i)
    {
        double momentum = kinematics.momenta[i-1];

        products[i]->SetDirection(GenerateRandomDirection());
        products[i]->SetMomentum(momentum);

        // Boost previous particles to new frame

        double M    = kinematics.virtual_masses[i - 1];
        double beta = momentum / std::sqrt(momentum * momentum + M * M);

        for (unsigned int s = 0; s < i; ++s)
        {
            // Boost in -p_i direction
            Boost(*products[s], products[i]->GetDirection(), beta);
        }
    }

    // Boost all daughters to parent frame
    Boost(products, particle.GetDirection(), particle.GetMomentum() / particle.GetEnergy());

    CopyParticleProperties(products, particle);

    return products;
}

// ------------------------------------------------------------------------- //
ManyBodyPhaseSpace::PhaseSpaceParameters ManyBodyPhaseSpace::GetPhaseSpaceParams(const ParticleDef& parent)
{
    ParameterMap::iterator it = parameter_map_.find(parent);

    if (it != parameter_map_.end())
    {
        return it->second;
    } else
    {
        PhaseSpaceParameters params;
        params.normalization = CalculateNormalization(parent.mass);
        params.weight_max = EstimateMaxWeight(params.normalization, parent.mass);
        parameter_map_[parent] = params;

        return params;
    }
}

// ------------------------------------------------------------------------- //
double ManyBodyPhaseSpace::CalculateNormalization(double parent_mass)
{
    double normalization = std::pow(parent_mass - sum_daughter_masses_, number_of_daughters_ - 2) /
                           boost::math::factorial<double>(number_of_daughters_ - 2) *
                           std::pow(2.0 * boost::math::constants::pi<double>(), number_of_daughters_ - 1.0);

    return normalization /= 2.0 * parent_mass;
}

// ------------------------------------------------------------------------- //
double ManyBodyPhaseSpace::EstimateMaxWeight(double normalization, double parent_mass)
{
    double weight = 1.0;
    double E_max = parent_mass - sum_daughter_masses_ + daughter_masses_[0];
    double E_min = 0.0;

    for (int i = 1; i < number_of_daughters_; ++i)
    {
        E_min += daughter_masses_[i - 1];
        E_max += daughter_masses_[i];
        weight *= Momentum(E_max, E_min, daughter_masses_[i]);
    }

    return normalization * weight;
}

// ------------------------------------------------------------------------- //
ManyBodyPhaseSpace::PhaseSpaceKinematics ManyBodyPhaseSpace::CalculateKinematics(double normalization, double parent_mass)
{
    PhaseSpaceKinematics kinematics;

    // Create sorted random numbers
    std::vector<double> randoms;
    randoms.reserve(daughters_.size());

    randoms.push_back(0.0);

    for (unsigned int i = 0; i < daughter_masses_.size() - 2; ++i)
    {
        randoms.push_back(RandomGenerator::Get().RandomDouble());
    }

    randoms.push_back(1.0);

    std::sort(randoms.begin(), randoms.end());

    // Calculate virtual masses
    double intermediate_mass = 0.0;
    for (int i = 0; i < number_of_daughters_; ++i)
    {
        intermediate_mass += daughter_masses_[i];
        kinematics.virtual_masses.push_back(intermediate_mass + randoms[i] * (parent_mass - sum_daughter_masses_));
    }

    // Calculate intermediate momenta
    double weight = 1.0;
    double momentum = 0.0;

    for (int i = 1; i < number_of_daughters_; ++i)
    {
        momentum = Momentum(kinematics.virtual_masses[i], kinematics.virtual_masses[i - 1], daughter_masses_[i]);
        kinematics.momenta.push_back(momentum);
        weight *= momentum;
    }

    kinematics.weight = normalization * weight;

    return kinematics;
}

// ------------------------------------------------------------------------- //
// Print
// ------------------------------------------------------------------------- //

void ManyBodyPhaseSpace::print(std::ostream& os) const
{
    os << "Used particle definitions:" << '\n';
    for (std::vector<const ParticleDef*>::const_iterator it = daughters_.begin(); it != daughters_.end(); ++it)
    {
        os << (*it)->name << '\t' << (*it)->mass << '\n';
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
    : daughters_(builder.daughters_)
{
}

ManyBodyPhaseSpace::Builder::~Builder()
{
    for (std::vector<const ParticleDef*>::const_iterator it = daughters_.begin(); it != daughters_.end(); ++it)
    {
        delete *it;
    }

    daughters_.clear();
}

// ------------------------------------------------------------------------- //
ManyBodyPhaseSpace::Builder& ManyBodyPhaseSpace::Builder::addDaughter(const ParticleDef& daughter)
{
    daughters_.push_back(daughter.clone());
    return *this;
}

// ------------------------------------------------------------------------- //
ManyBodyPhaseSpace ManyBodyPhaseSpace::Builder::build()
{
    return ManyBodyPhaseSpace(daughters_);
}
