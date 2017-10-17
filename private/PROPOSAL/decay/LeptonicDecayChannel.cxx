
#include <boost/bind.hpp>
#include <cmath>

#include <boost/math/tools/roots.hpp>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/decay/LeptonicDecayChannel.h"
#include "PROPOSAL/particle/Particle.h"

#include "PROPOSAL/Output.h"

using namespace PROPOSAL;

const std::string LeptonicDecayChannel::name_ = "LeptonicDecayChannel";

LeptonicDecayChannel::LeptonicDecayChannel(const ParticleDef& lepton,
                                           const ParticleDef& neutrino,
                                           const ParticleDef& anti_neutrino)
    : DecayChannel()
    , massive_lepton_(lepton)
    // , neutrino_()
    // , anti_neutrino_()
    , neutrino_(neutrino)
    , anti_neutrino_(anti_neutrino)
    // , neutrino_(GetNeutrino(massive_lepton_))
    // , anti_neutrino_(GetAntiNeutrino(massive_lepton_))
    , root_finder_(IMAXS, IPREC)
{
}

LeptonicDecayChannel::~LeptonicDecayChannel()
{
}

LeptonicDecayChannel::LeptonicDecayChannel(const LeptonicDecayChannel& mode)
    : DecayChannel(mode)
    , massive_lepton_(mode.massive_lepton_)
    , neutrino_(mode.neutrino_)
    , anti_neutrino_(mode.anti_neutrino_)
    , root_finder_(mode.root_finder_)
{
}

bool LeptonicDecayChannel::compare(const DecayChannel& channel) const
{
    const LeptonicDecayChannel* leptonic = dynamic_cast<const LeptonicDecayChannel*>(&channel);

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

double LeptonicDecayChannel::DecayRate(double x, double right_side)
{
    double x2;
    x2  =   x*x;
    return  x*x2 - x2*x2/2 - right_side;
}

double LeptonicDecayChannel::DifferentialDecayRate(double x)
{
    return (3 - 2*x) * x*x;
}

std::pair<double, double> LeptonicDecayChannel::function_and_derivative(double x, double right_side)
{
    return std::make_pair (DecayRate(x, right_side), DifferentialDecayRate(x));
}

double LeptonicDecayChannel::FindRootBoost(double min, double right_side)
{
    double max = 1;
    double x_start = 0.5;
    int binary_digits = 6;
    // int max_steps = 20;

    return boost::math::tools::newton_raphson_iterate(
        boost::bind(&LeptonicDecayChannel::function_and_derivative, this, _1, right_side)
        , x_start
        , min
        , max
        , binary_digits
        );
}


DecayChannel::DecayProducts LeptonicDecayChannel::Decay(Particle& particle)
{
    double emax, x0, f0, el, pl, right_side;
    double lm  =   ME;
    double parent_mass = particle.GetMass();

    DecayProducts products;
    products.push_back(new Particle(massive_lepton_));
    products.push_back(new Particle(neutrino_));
    products.push_back(new Particle(anti_neutrino_));

    // Sample energy from decay rate
    double ernd = RandomGenerator::Get().RandomDouble();
    // double arnd = RandomGenerator::Get().RandomDouble();

    emax    =   (parent_mass*parent_mass + lm*lm) / (2*parent_mass);
    x0      =   lm/emax;
    f0      =   x0*x0;
    f0      =   f0*x0 - f0*f0/2;
    right_side = f0+(0.5-f0)*ernd;

    // double find_root_old = root_finder_.FindRoot(
    //         x0, 1.0, 0.5,
    //         boost::bind(&LeptonicDecayChannel::DecayRate, this, _1, right_side),
    //         boost::bind(&LeptonicDecayChannel::DifferentialDecayRate, this, _1));

    double find_root_new = FindRootBoost(x0, right_side);

    // printf("%f %f\n", find_root_old, find_root_new);

    el = std::max(find_root_new *emax, lm);
    pl = sqrt(el*el - lm*lm);

    // Sample directions For the massive letpon
    products[0]->SetDirection(GenerateRandomDirection());
    products[0]->SetMomentum(pl);

    // Sample directions For the massless letpon
    double energy_neutrinos = parent_mass - el;
    double virtual_mass     = std::sqrt((energy_neutrinos - pl) * (energy_neutrinos + pl));
    double momentum_neutrinos = 0.5 * virtual_mass;
    Vector3D direction = GenerateRandomDirection();

    products[1]->SetDirection(direction);
    products[1]->SetMomentum(momentum_neutrinos);

    products[2]->SetDirection(-direction);
    products[2]->SetMomentum(momentum_neutrinos);

    // Boost neutrinos to lepton frame
    double beta = pl / energy_neutrinos;
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

void LeptonicDecayChannel::print(std::ostream& os) const
{
    os << "Massive lepton:\n" << massive_lepton_ << '\n';
    os << "Neutrino:\n" << neutrino_ << '\n';
    os << "Anti neutrino:\n" << anti_neutrino_ << '\n';
}
