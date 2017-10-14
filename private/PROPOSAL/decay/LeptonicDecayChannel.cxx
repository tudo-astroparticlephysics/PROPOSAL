
#include <boost/bind.hpp>
#include <cmath>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/MathModel.h"
#include "PROPOSAL/decay/LeptonicDecayChannel.h"
#include "PROPOSAL/particle/PROPOSALParticle.h"

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

double LeptonicDecayChannel::DecayRate(double x)
{
    double x2;
    x2  =   x*x;
    return  x*x2 - x2*x2/2;
}

double LeptonicDecayChannel::DifferentialDecayRate(double x)
{
    return (3 - 2*x) * x*x;
}

DecayChannel::DecayProducts LeptonicDecayChannel::Decay(PROPOSALParticle* particle)
{
    double emax, x0, f0, el, pl, final_energy;
    double lm  =   ME;
    double parent_mass = particle->GetMass();

    double ernd = RandomGenerator::Get().RandomDouble();
    double arnd = RandomGenerator::Get().RandomDouble();

    emax    =   (parent_mass*parent_mass + lm*lm) / (2*parent_mass);
    x0      =   lm/emax;
    f0      =   x0*x0;
    f0      =   f0*x0 - f0*f0/2;
    el      =   std::max(root_finder_.FindRoot(
                        x0, 1.0, 0.5,
                        // &LeptonicDecayChannel::DecayRate,
                        // &LeptonicDecayChannel::DifferentialDecayRate,
                        boost::bind(&LeptonicDecayChannel::DecayRate, this, _1),
                        boost::bind(&LeptonicDecayChannel::DifferentialDecayRate, this, _1),
                        f0+(0.5-f0)*ernd) *emax
                    , lm);
    pl      =   sqrt(el*el - lm*lm);

    final_energy = el * (particle->GetEnergy()/parent_mass) + pl * (particle->GetMomentum()/parent_mass) * (2*arnd - 1);

    // Create products
    PROPOSALParticle* product_particle = new PROPOSALParticle(ParticleDef::Builder().SetEMinus().build());
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
