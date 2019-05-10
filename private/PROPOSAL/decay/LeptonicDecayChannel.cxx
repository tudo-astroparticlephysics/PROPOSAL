
#include <functional>
#include <cmath>

#include <boost/math/tools/roots.hpp>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/decay/LeptonicDecayChannel.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/particle/Particle.h"


using namespace PROPOSAL;

/******************************************************************************
*                         LeptonicDecayChannelApprox                         *
******************************************************************************/


const std::string LeptonicDecayChannelApprox::name_ = "LeptonicDecayChannelApprox";

// ------------------------------------------------------------------------- //
LeptonicDecayChannelApprox::LeptonicDecayChannelApprox(const ParticleDef& lepton,
                                           const ParticleDef& neutrino,
                                           const ParticleDef& anti_neutrino)
    : DecayChannel()
    , massive_lepton_(lepton)
    , neutrino_(neutrino)
    , anti_neutrino_(anti_neutrino)
{
}

// ------------------------------------------------------------------------- //
LeptonicDecayChannelApprox::~LeptonicDecayChannelApprox() {}

// ------------------------------------------------------------------------- //
LeptonicDecayChannelApprox::LeptonicDecayChannelApprox(const LeptonicDecayChannelApprox& mode)
    : DecayChannel(mode)
    , massive_lepton_(mode.massive_lepton_)
    , neutrino_(mode.neutrino_)
    , anti_neutrino_(mode.anti_neutrino_)
{
}

// ------------------------------------------------------------------------- //
bool LeptonicDecayChannelApprox::compare(const DecayChannel& channel) const
{
    const LeptonicDecayChannelApprox* leptonic = dynamic_cast<const LeptonicDecayChannelApprox*>(&channel);

    if (!leptonic)
        return false;
    else if (massive_lepton_ != leptonic->massive_lepton_)
        return false;
    else if (neutrino_ != leptonic->neutrino_)
        return false;
    else if (anti_neutrino_ != leptonic->anti_neutrino_)
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
double LeptonicDecayChannelApprox::DecayRate(double x, double parent_mass, double E_max, double right_side)
{
    (void)parent_mass;
    (void)E_max;

    return x * x * x * (1. - 0.5 * x) - right_side;
}

// ------------------------------------------------------------------------- //
double LeptonicDecayChannelApprox::DifferentialDecayRate(double x, double parent_mass, double E_max)
{
    (void)parent_mass;
    (void)E_max;

    return (3 - 2 * x) * x * x;
}

// ------------------------------------------------------------------------- //
std::pair<double, double> LeptonicDecayChannelApprox::function_and_derivative(double x,
                                                                        double parent_mass,
                                                                        double E_max,
                                                                        double right_side)
{
    return std::make_pair(DecayRate(x, parent_mass, E_max, right_side), DifferentialDecayRate(x, parent_mass, E_max));
}

// ------------------------------------------------------------------------- //
double LeptonicDecayChannelApprox::FindRootBoost(double min, double parent_mass, double E_max, double right_side)
{
    double max        = 1;
    double x_start    = 0.5;
    int binary_digits = 6;
    // in older versions a max_step was set to 40, which were the max number of int steps
    // int max_steps = 40;

    return boost::math::tools::newton_raphson_iterate(
        std::bind(&LeptonicDecayChannelApprox::function_and_derivative, this, std::placeholders::_1, parent_mass, E_max, right_side),
        x_start,
        min,
        max,
        binary_digits);
}

// ------------------------------------------------------------------------- //
DecayChannel::DecayProducts LeptonicDecayChannelApprox::Decay(const Particle& particle)
{
    double parent_mass = particle.GetMass();

    DecayProducts products;
    products.push_back(new Particle(massive_lepton_));
    products.push_back(new Particle(neutrino_));
    products.push_back(new Particle(anti_neutrino_));

    // Sample energy from decay rate
    double emax       = (parent_mass * parent_mass + massive_lepton_.mass * massive_lepton_.mass) / (2 * parent_mass);
    double x_min      = massive_lepton_.mass / emax;
    // double f_min      = x_min * x_min * x_min * (1 - 0.5 * x_min);
    // double right_side = f_min + (0.5 - f_min) * RandomGenerator::Get().RandomDouble();

    double f_min      = DecayRate(x_min, parent_mass, emax, 0.0);
    double f_max      = DecayRate(1.0, parent_mass, emax, 0.0);
    double right_side = f_min + (f_max - f_min) * RandomGenerator::Get().RandomDouble();

    double find_root = FindRootBoost(x_min, parent_mass, emax, right_side);

    double lepton_energy   = std::max(find_root * emax, massive_lepton_.mass);
    double lepton_momentum = sqrt(lepton_energy * lepton_energy - massive_lepton_.mass * massive_lepton_.mass);

    // Sample directions For the massive letpon
    products[0]->SetDirection(GenerateRandomDirection());
    products[0]->SetMomentum(lepton_momentum);

    // Sample directions For the massless letpon
    double energy_neutrinos   = parent_mass - lepton_energy;
    double virtual_mass       = std::sqrt((energy_neutrinos - lepton_momentum) * (energy_neutrinos + lepton_momentum));
    double momentum_neutrinos = 0.5 * virtual_mass;
    Vector3D direction        = GenerateRandomDirection();

    products[1]->SetDirection(direction);
    products[1]->SetMomentum(momentum_neutrinos);

    products[2]->SetDirection(-direction);
    products[2]->SetMomentum(momentum_neutrinos);

    // Boost neutrinos to lepton frame
    double beta = lepton_momentum / energy_neutrinos;
    Boost(*products[1], products[0]->GetDirection(), -beta);
    Boost(*products[2], products[0]->GetDirection(), -beta);

    // Boost all particle in parent frame
    Boost(products, particle.GetDirection(), particle.GetMomentum() / particle.GetEnergy());

    CopyParticleProperties(products, particle);

    return products;
}

// ------------------------------------------------------------------------- //
// Print
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
void LeptonicDecayChannelApprox::print(std::ostream& os) const
{
    os << "Massive lepton:\n" << massive_lepton_ << '\n';
    os << "Neutrino:\n" << neutrino_ << '\n';
    os << "Anti neutrino:\n" << anti_neutrino_ << '\n';
}

/******************************************************************************
 *                          LeptonicDecayChannel                              *
 ******************************************************************************/

const std::string LeptonicDecayChannel::name_ = "LeptonicDecayChannel";

// ------------------------------------------------------------------------- //
LeptonicDecayChannel::LeptonicDecayChannel(const ParticleDef& lepton,
                                                 const ParticleDef& neutrino,
                                                 const ParticleDef& anti_neutrino)
    : LeptonicDecayChannelApprox(lepton, neutrino, anti_neutrino)
{
}

// ------------------------------------------------------------------------- //
LeptonicDecayChannel::~LeptonicDecayChannel() {}

// ------------------------------------------------------------------------- //
LeptonicDecayChannel::LeptonicDecayChannel(const LeptonicDecayChannel& mode)
    : LeptonicDecayChannelApprox(mode)
{
}

// ------------------------------------------------------------------------- //
double LeptonicDecayChannel::DecayRate(double x, double M, double E_max, double right_side)
{
    double M2 = M * M;
    double m  = massive_lepton_.mass;
    double m2 = m * m;

    double E_l     = E_max * x;
    double sqrt_EM = std::sqrt(E_l * E_l - m2);

    return 1.5 * m2 * m2 * M * std::log(sqrt_EM + E_l) +
           sqrt_EM * ((M2 + m2 - M * E_l) * (E_l * E_l - m2) - 1.5 * M * E_l * m2) - right_side;
}

// ------------------------------------------------------------------------- //
double LeptonicDecayChannel::DifferentialDecayRate(double x, double M, double E_max)
{
    double m   = massive_lepton_.mass;
    double E_l = E_max * x;

    return E_max * std::sqrt(E_l * E_l - m * m) * (M * E_l * (3.0 * M - 4.0 * E_l) + m * m * (3.0 * E_l - 2 * M));
}
