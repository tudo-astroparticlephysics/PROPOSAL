
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
ManyBodyPhaseSpace::ManyBodyPhaseSpace(std::vector<std::shared_ptr<const ParticleDef>> daughters, MatrixElementFunction me)
    : DecayChannel()
    , daughters_(daughters)
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
    for (const auto& i : daughters) {

        daughter_masses_.push_back(i->mass);
        sum_daughter_masses_ += i->mass;
    }
}


// ------------------------------------------------------------------------- //
ManyBodyPhaseSpace::ManyBodyPhaseSpace(const ManyBodyPhaseSpace& mode)
    : DecayChannel(mode)
    , daughters_(mode.daughters_)
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
std::vector<ParticleState> ManyBodyPhaseSpace::Decay(const ParticleDef& p_def, const ParticleState& p_condition)
{
    // Create vector for decay products
    std::vector<ParticleState> products;

    for (auto p : daughters_) {
        products.emplace_back((ParticleType)p->particle_type, p_condition.position, p_condition.direction, p_condition.energy, p_condition.time, 0);
    }

    // prefactor for the phase space density
    PhaseSpaceParameters params = GetPhaseSpaceParams(p_def);
    PhaseSpaceKinematics kinematics;

    if (uniform_)
    {
        double weight_ref;
        double weight_sample;
        do
        {
            // precalculated kinematics
            kinematics = CalculateKinematics(params.normalization, p_def.mass);
            GenerateEvent(products, kinematics);
            // sample product states with rejection sampling
            weight_ref = params.weight_min + RandomGenerator::Get().RandomDouble() * (params.weight_max - params.weight_min);
            weight_sample = kinematics.weight * matrix_element_(p_condition, products);

        } while(weight_ref > weight_sample);
    }
    else
    {
        // precalculated kinematics
        kinematics = CalculateKinematics(params.normalization, p_def.mass);
        GenerateEvent(products, kinematics);
    }

    // Boost all products in Lab frame (the reason, why the boosting goes in the negative direction of the particle)
    Boost(products, -p_condition.direction, p_condition.energy/p_def.mass, p_condition.GetMomentum() / p_def.mass);

    return products;
}

// ------------------------------------------------------------------------- //
void ManyBodyPhaseSpace::GenerateEvent(std::vector<ParticleState>& products, const PhaseSpaceKinematics& kinematics)
{
    // Calculate first momentum in R2
    Cartesian3D direction = GenerateRandomDirection();

    products[1].direction = direction;
    products[1].SetMomentum(kinematics.momenta[0]);

    Cartesian3D opposite_direction = -direction;
    products[0].direction = opposite_direction;
    products[0].SetMomentum(kinematics.momenta[0]);

    // Correct the previous momenta
    for (unsigned int i = 2; i < daughter_masses_.size(); ++i)
    {
        double momentum = kinematics.momenta[i-1];

        products[i].direction = GenerateRandomDirection();
        products[i].SetMomentum(momentum);

        // Boost previous particles to new frame

        double M = kinematics.virtual_masses[i - 1];
        // double beta = momentum / std::sqrt(momentum * momentum + M * M);
        double gamma = std::sqrt(momentum * momentum + M * M) / M;
        double betagamma = momentum / M;

        for (unsigned int s = 0; s < i; ++s)
        {
            // Boost in -p_i direction
            Boost(products[s], products[i].direction, gamma, betagamma);
        }
    }
}

// ------------------------------------------------------------------------- //
double ManyBodyPhaseSpace::DefaultEvaluate(const ParticleState& p_condition, const std::vector<ParticleState>& products)
{
    (void) p_condition;
    (void) products;

    return 1.0;
}

// ------------------------------------------------------------------------- //
double ManyBodyPhaseSpace::Evaluate(const ParticleState& p_condition, const std::vector<ParticleState>& products)
{
    return matrix_element_(p_condition, products);
}

// ------------------------------------------------------------------------- //
ManyBodyPhaseSpace::PhaseSpaceParameters ManyBodyPhaseSpace::GetPhaseSpaceParams(const ParticleDef& parent_def)
{
    ParameterMap::iterator it = parameter_map_.find(parent_def);

    if (it != parameter_map_.end())
    {
        return it->second;
    } else
    {
        PhaseSpaceParameters params;

        params.normalization = CalculateNormalization(parent_def.mass);
        estimate_(params, parent_def);

        parameter_map_[parent_def] = params;

        return params;
    }
}

// ------------------------------------------------------------------------- //
double ManyBodyPhaseSpace::CalculateNormalization(double parent_mass)
{
    // For mario: n! = Gamma(n + 1)
    double normalization = std::pow(parent_mass - sum_daughter_masses_, number_of_daughters_ - 2) /
                           std::tgamma(number_of_daughters_ - 2.0 + 1.0) *
                           std::pow(2.0 * PI, number_of_daughters_ - 1.0);

    return normalization / (2.0 * parent_mass);
}

// ------------------------------------------------------------------------- //
void ManyBodyPhaseSpace::EstimateMaxWeight(PhaseSpaceParameters& params, const ParticleDef& parent_def)
{
    double weight = 1.0;
    double E_max = parent_def.mass - sum_daughter_masses_ + daughter_masses_[0];
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
void ManyBodyPhaseSpace::SampleEstimateMaxWeight(PhaseSpaceParameters& params, const ParticleDef& parent_def)
{
    // Create vector for decay products
    std::vector<ParticleState> products;

    for (auto d : daughters_) {
        products.emplace_back();
        products.back().type = d->particle_type;
    }

    // precalculated kinematics
    ParticleState particle;
    particle.type = parent_def.particle_type;
    particle.energy = parent_def.mass;

    // initialization of weights
    PhaseSpaceKinematics kinematics = CalculateKinematics(params.normalization, parent_def.mass);
    GenerateEvent(products, kinematics);
    double result = kinematics.weight * matrix_element_(particle, products);
    params.weight_min = result;
    params.weight_max = result;

    for (int i = 1; i < broad_phase_statistic_; ++i)
    {
        kinematics = CalculateKinematics(params.normalization, parent_def.mass);
        GenerateEvent(products, kinematics);
        result = kinematics.weight * matrix_element_(particle, products);

        if (result < params.weight_min)
        {
            params.weight_min = result; // the lower limit is at the edge of the phase space, which is zero, or?
        }

        if (result > params.weight_max)
        {
            params.weight_max = result;
        }
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
    for (const auto& it : daughters_)
        os << it->name << '\t' << it->mass << '\n';
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


// ------------------------------------------------------------------------- //
ManyBodyPhaseSpace::Builder& ManyBodyPhaseSpace::Builder::addDaughter(const ParticleDef& daughter)
{
    daughters_.push_back(std::make_shared<ParticleDef>(daughter));
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
