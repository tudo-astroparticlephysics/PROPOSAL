
// #include "PROPOSAL/Constants.h"

#include <algorithm> // std::sort
#include <cmath>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/decay/ManyBodyPhaseSpace.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/particle/Particle.h"


using namespace PROPOSAL;

const std::string ManyBodyPhaseSpace::name_ = "ManyBodyPhaseSpace";

// ------------------------------------------------------------------------- //
ManyBodyPhaseSpace::ManyBodyPhaseSpace(std::vector<const ParticleDef*> daughters, MatrixElementFunction me)
    : DecayChannel()
    , daughters_()
    , daughter_masses_()
    , number_of_daughters_(static_cast<int>(daughters.size()))
    , sum_daughter_masses_(0.0)
    , uniform_(true)
    , broad_phase_statistic_(1000)
    , matrix_element_()
    , use_default_matrix_element_(true)
    , estimate_(nullptr)
{
    if (me == nullptr)
    {
        matrix_element_ = ManyBodyPhaseSpace::DefaultEvaluate;
        use_default_matrix_element_ = true;
        estimate_ = std::bind(&ManyBodyPhaseSpace::EstimateMaxWeight, this, std::placeholders::_1, std::placeholders::_2);
    }
    else
    {
        matrix_element_ = me;
        use_default_matrix_element_ = false;
        estimate_ = std::bind(&ManyBodyPhaseSpace::SampleEstimateMaxWeight, this, std::placeholders::_1, std::placeholders::_2);
    }

    for (std::vector<const ParticleDef*>::iterator it = daughters.begin(); it != daughters.end(); ++it)
    {
        daughters_.push_back((*it)->clone());
        daughter_masses_.push_back((*it)->mass);
        sum_daughter_masses_ += (*it)->mass;
    }
}

// ------------------------------------------------------------------------- //
ManyBodyPhaseSpace::~ManyBodyPhaseSpace()
{
    for (std::vector<const ParticleDef*>::const_iterator it = daughters_.begin(); it != daughters_.end(); ++it)
    {
        delete *it;
    }

    daughters_.clear();
}

// ------------------------------------------------------------------------- //
ManyBodyPhaseSpace::ManyBodyPhaseSpace(const ManyBodyPhaseSpace& mode)
    : DecayChannel(mode)
    , daughters_()
    , daughter_masses_(mode.daughter_masses_)
    , number_of_daughters_(mode.number_of_daughters_)
    , sum_daughter_masses_(mode.sum_daughter_masses_)
    , uniform_(mode.uniform_)
    , broad_phase_statistic_(mode.broad_phase_statistic_)
    , matrix_element_(mode.matrix_element_)
    , use_default_matrix_element_(mode.use_default_matrix_element_)
{
    if (use_default_matrix_element_)
    {
        estimate_ = std::bind(&ManyBodyPhaseSpace::EstimateMaxWeight, this, std::placeholders::_1, std::placeholders::_2);
    }
    else
    {
        estimate_ = std::bind(&ManyBodyPhaseSpace::SampleEstimateMaxWeight, this, std::placeholders::_1, std::placeholders::_2);
    }

    for (std::vector<const ParticleDef*>::const_iterator it = mode.daughters_.begin(); it != mode.daughters_.end();
         ++it)
    {
        daughters_.push_back((*it)->clone());
    }
}

// ------------------------------------------------------------------------- //
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
    } else if (broad_phase_statistic_ != many_body->broad_phase_statistic_)
    {
        return false;
    } else if (use_default_matrix_element_ != many_body->use_default_matrix_element_)
    {
        return false;
    } else
    {
        return true;
    }
}

// ------------------------------------------------------------------------- //
DecayChannel::DecayProducts ManyBodyPhaseSpace::Decay(const Particle& particle)
{
    // Create vector for decay products
    DecayProducts products;
    products.reserve(daughters_.size());

    for (std::vector<const ParticleDef*>::const_iterator it = daughters_.begin(); it != daughters_.end(); ++it)
    {
        products.push_back(new Particle(**it));
    }

    // prefactor for the phase space density
    PhaseSpaceParameters params = GetPhaseSpaceParams(particle.GetParticleDef());

    if (uniform_)
    {
        double weight = 0.0;

        do
        {
            weight = GenerateEvent(products, params, particle);

        } while(params.weight_min + RandomGenerator::Get().RandomDouble() * (params.weight_max - params.weight_min) > weight * matrix_element_(particle, products));
    }
    else
    {
        GenerateEvent(products, params, particle);
    }

    // Boost all products in Lab frame (the reason, why the boosting goes in the negative direction of the particle)
    Boost(products, -particle.GetDirection(), particle.GetEnergy()/particle.GetMass(), particle.GetMomentum() / particle.GetMass());

    CopyParticleProperties(products, particle);

    return products;
}

// ------------------------------------------------------------------------- //
double ManyBodyPhaseSpace::GenerateEvent(DecayProducts& products, const PhaseSpaceParameters& params, const Particle& particle)
{
    double parent_mass = particle.GetMass();

    // precalculated kinematics
    PhaseSpaceKinematics kinematics = CalculateKinematics(params.normalization, parent_mass);

    // Calculate first momentum in R2
    Vector3D direction = GenerateRandomDirection();

    products[1]->SetDirection(direction);
    products[1]->SetMomentum(kinematics.momenta[0]);

    Vector3D opposite_direction = -direction;
    opposite_direction.CalculateSphericalCoordinates();
    products[0]->SetDirection(opposite_direction);
    products[0]->SetMomentum(kinematics.momenta[0]);

    // Correct the previous momenta
    for (unsigned int i = 2; i < daughter_masses_.size(); ++i)
    {
        double momentum = kinematics.momenta[i-1];

        products[i]->SetDirection(GenerateRandomDirection());
        products[i]->SetMomentum(momentum);

        // Boost previous particles to new frame

        double M    = kinematics.virtual_masses[i - 1];
        // double beta = momentum / std::sqrt(momentum * momentum + M * M);
        double gamma = std::sqrt(momentum * momentum + M * M) / M;
        double betagamma = momentum / M;

        for (unsigned int s = 0; s < i; ++s)
        {
            // Boost in -p_i direction
            Boost(*products[s], products[i]->GetDirection(), gamma, betagamma);
        }
    }

    return kinematics.weight;
}

// ------------------------------------------------------------------------- //
double ManyBodyPhaseSpace::DefaultEvaluate(const Particle& particle, const ManyBodyPhaseSpace::DecayProducts& products)
{
    (void) particle;
    (void) products;

    return 1.0;
}

// ------------------------------------------------------------------------- //
double ManyBodyPhaseSpace::Evaluate(const Particle& particle, const ManyBodyPhaseSpace::DecayProducts& products)
{
    return matrix_element_(particle, products);
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

        CalculateNormalization(params, parent.mass);
        estimate_(params, parent);

        parameter_map_[parent] = params;

        return params;
    }
}

// ------------------------------------------------------------------------- //
void ManyBodyPhaseSpace::CalculateNormalization(PhaseSpaceParameters& params, double parent_mass)
{
    // For mario: n! = Gamma(n + 1)
    double normalization = std::pow(parent_mass - sum_daughter_masses_, number_of_daughters_ - 2) /
                           std::tgamma(number_of_daughters_ - 2.0 + 1.0) *
                           std::pow(2.0 * PI, number_of_daughters_ - 1.0);

    params.normalization = normalization /= 2.0 * parent_mass;
}

// ------------------------------------------------------------------------- //
void ManyBodyPhaseSpace::EstimateMaxWeight(PhaseSpaceParameters& params, const ParticleDef& parent)
{
    double weight = 1.0;
    double E_max = parent.mass - sum_daughter_masses_ + daughter_masses_[0];
    double E_min = 0.0;

    for (int i = 1; i < number_of_daughters_; ++i)
    {
        E_min += daughter_masses_[i - 1];
        E_max += daughter_masses_[i];
        weight *= Momentum(E_max, E_min, daughter_masses_[i]);
    }

    params.weight_max = params.normalization * weight;
    params.weight_min = 0.0;
}

// ------------------------------------------------------------------------- //
void ManyBodyPhaseSpace::SampleEstimateMaxWeight(PhaseSpaceParameters& params, const ParticleDef& parent)
{
    // Create vector for decay products
    DecayProducts products;
    products.reserve(daughters_.size());

    for (std::vector<const ParticleDef*>::const_iterator it = daughters_.begin(); it != daughters_.end(); ++it)
    {
        products.push_back(new Particle(**it));
    }

    Particle particle(parent);

    // precalculated kinematics
    PhaseSpaceKinematics kinematics = CalculateKinematics(params.normalization, parent.mass);

    double result = 0.0;

    for (int i = 0; i < broad_phase_statistic_; ++i)
    {
        result = GenerateEvent(products, params, particle) * matrix_element_(particle, products);

        if (result < params.weight_min)
        {
            params.weight_min = result;
        }

        if (result > params.weight_max)
        {
            params.weight_max = result;
        }
    }

    for (DecayProducts::const_iterator it = products.begin(); it != products.end(); ++it)
    {
        delete *it;
    }
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
    , matrix_element_()
{
}

ManyBodyPhaseSpace::Builder::Builder(const Builder& builder)
    : daughters_(builder.daughters_)
    , matrix_element_(builder.matrix_element_)
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
ManyBodyPhaseSpace::Builder& ManyBodyPhaseSpace::Builder::setMatrixElement(MatrixElementFunction me)
{
    matrix_element_ = me;
    return *this;
}

// ------------------------------------------------------------------------- //
ManyBodyPhaseSpace ManyBodyPhaseSpace::Builder::build()
{
    return ManyBodyPhaseSpace(daughters_, matrix_element_);
}
