
#include <boost/bind.hpp>
#include <cmath>

#include <boost/math/tools/roots.hpp>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/decay/LeptonicDecayChannel.h"
#include "PROPOSAL/particle/Particle.h"

using namespace PROPOSAL;

const std::string LeptonicDecayChannel::name_ = "LeptonicDecayChannel";

LeptonicDecayChannel::LeptonicDecayChannel()
    : DecayChannel()
    , root_finder_(IMAXS, IPREC)
{
}

LeptonicDecayChannel::~LeptonicDecayChannel()
{
}

LeptonicDecayChannel::LeptonicDecayChannel(const LeptonicDecayChannel& mode)
    : DecayChannel(mode)
    , root_finder_(mode.root_finder_)
{
}

bool LeptonicDecayChannel::compare(const DecayChannel& channel) const
{
    const LeptonicDecayChannel* leptonic = dynamic_cast<const LeptonicDecayChannel*>(&channel);

    if (!leptonic)
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


DecayChannel::DecayProducts LeptonicDecayChannel::Decay(Particle* particle)
{
    double emax, x0, f0, el, pl, final_energy, right_side;
    double lm  =   ME;
    double parent_mass = particle->GetMass();

    double ernd = RandomGenerator::Get().RandomDouble();
    double arnd = RandomGenerator::Get().RandomDouble();

    emax    =   (parent_mass*parent_mass + lm*lm) / (2*parent_mass);
    x0      =   lm/emax;
    f0      =   x0*x0;
    f0      =   f0*x0 - f0*f0/2;
    right_side = f0+(0.5-f0)*ernd;

    double find_root_old = root_finder_.FindRoot(
            x0, 1.0, 0.5,
            boost::bind(&LeptonicDecayChannel::DecayRate, this, _1, right_side),
            boost::bind(&LeptonicDecayChannel::DifferentialDecayRate, this, _1));

    double find_root_new = FindRootBoost(x0, right_side);

    printf("%f %f\n", find_root_old, find_root_new);

    el = std::max(find_root_new *emax, lm);
    pl = sqrt(el*el - lm*lm);

    final_energy = el * (particle->GetEnergy()/parent_mass) + pl * (particle->GetMomentum()/parent_mass) * (2*arnd - 1);

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

void LeptonicDecayChannel::print(std::ostream& os) const
{
    os << "Lepton mass: " << '\n';
}
