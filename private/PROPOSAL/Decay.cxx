/*
* Decay.cxx
*
*  Created on: 24.06.2010
*      Author: koehne
*/

// #include <algorithm>
// #include <cmath>

#include <boost/bind.hpp>

#include "PROPOSAL/Decay.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/MathModel.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Output.h"

using namespace std;
using namespace PROPOSAL;


/******************************************************************************
*                                   Decay                                     *
******************************************************************************/

double Decay::MakeDecay(PROPOSALParticle* particle)
{

    if(multiplier_ <= 0 || particle->GetLifetime() < 0)
    {
        return 0;
    }

    return multiplier_/max((particle->GetMomentum()/particle->GetMass())*particle->GetLifetime()*SPEED, XRES);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// double Decay::CalculateProductEnergy( double ernd, double arnd, double srnd )
// {
//     if(particle_->GetLifetime()<0)
//     {
//         // out_ = ParticleType::unknown;
//         return 0;
//     }
//     double emax, x0, f0, el, lm, pl;
//     ParticleType::Enum out1 = ParticleType::NuMu;
//     ParticleType::Enum out2 = ParticleType::NuMu;
//
//     // Tau Decay
//     if (particle_->GetType() == ParticleType::TauMinus || particle_->GetType() == ParticleType::TauPlus)
//     {
//
//         const double brmu   =   0.1737;         // ratio: tau- ---> mu + anti-munu + nu-tau
//         const double brel   =   0.1783+brmu;    // ratio: tau- ---> electron + anti-electron +nu-tau
//         const double brpi   =   0.1109+brel;    // ratio: tau- ---> pi- + nu-tau
//         const double br2p   =   0.2540+brpi;    // ratio: tau- ---> pi- + pi0 + nu-tau
//         const double br3p   =   0.1826+br2p;    // ratio: tau- ---> pi- + pi0 + pi0 + nutau and pi- + pi+ + pi- + nu-tau
//
//         if( srnd<brmu )
//         {
//             lm  =   MMU;
//
//             if (particle_->GetType() == ParticleType::TauPlus)
//             {
//                 out_ = ParticleType::MuPlus;
//                 if(store_neutrinos_)
//                 {
//                     out1 = ParticleType::NuTauBar;
//                     out2 = ParticleType::NuMu;
//                 }
//             }
//             else if (particle_->GetType() == ParticleType::TauMinus)
//             {
//                 out_ = ParticleType::MuMinus;
//                 if(store_neutrinos_)
//                 {
//                     out1 = ParticleType::NuTau;
//                     out2 = ParticleType::NuMuBar;
//                 }
//             }
//             else
//             {
//                 log_fatal("You should never end here: tau- or tau+?");
//             }
//         }
//         else if( srnd<brel )
//         {
//             lm  =   ME;
//
//             if (particle_->GetType() == ParticleType::TauPlus)
//             {
//                 out_ = ParticleType::EPlus;
//
//                 if(store_neutrinos_)
//                 {
//                     out1 = ParticleType::NuTauBar;
//                     out2 = ParticleType::NuE;
//                 }
//             }
//             else if (particle_->GetType() == ParticleType::TauMinus)
//             {
//                 out_ = ParticleType::EMinus;
//
//                 if(store_neutrinos_)
//                 {
//                     out1 = ParticleType::NuTau;
//                     out2 = ParticleType::NuEBar;
//                 }
//             }
//             else
//             {
//                 log_fatal("You should never end here: tau- or tau+?");
//             }
//         }
//         else
//         {
//             if( srnd<brpi )
//             {
//                 lm  =   MPI;
//             }
//             else if( srnd<br2p )
//             {
//                 lm  =   MRH;
//             }
//             else if( srnd<br3p )
//             {
//                 lm  =   MA1;
//             }
//             else
//             {
//                 lm  =   MRS;
//             }
//             el  =   (MTAU*MTAU + lm*lm) / (2*MTAU);
//
//             out_ = ParticleType::Hadrons;
//
//             el  =   el * (particle_->GetEnergy()/particle_->GetMass()) + sqrt(el*el - lm*lm) * (particle_->GetMomentum()/particle_->GetMass()) * (2*arnd-1);
//
//             if(store_neutrinos_)
//             {
//                 if( EndsWith(particle_->GetName(),"+") )
//                 {
//                     out1 = ParticleType::NuTauBar;
//                 }
//                 else if( EndsWith(particle_->GetName(),"-") )
//                 {
//                     out1 = ParticleType::NuTau;
//                 }
//                 else
//                 {
//                     out1 = ParticleType::NuTau;
//                 }
//                // o->output(1, out1, particle_->GetEnergy()-el, 0);
//             }
//             return el;
//         }
//     }
//
//     //Muon Decay
//     else
//     {
//         lm  =   ME;
//         if (particle_->GetType() == ParticleType::MuPlus)
//         {
//             out_ = ParticleType::EPlus;
//
//             if(store_neutrinos_)
//             {
//                 out1 = ParticleType::NuMuBar;
//                 out2 = ParticleType::NuE;
//             }
//         }
//         else if (particle_->GetType() == ParticleType::MuMinus)
//         {
//             out_ = ParticleType::EMinus;
//
//             if(store_neutrinos_)
//             {
//                 out1 = ParticleType::NuMu;
//                 out2 = ParticleType::NuEBar;
//             }
//         }
//         else
//         {
//             log_fatal("You should never end here: mu- or mu+?");
//         }
//     }
//     emax    =   (particle_->GetMass()*particle_->GetMass() + lm*lm) / (2*particle_->GetMass());
//     x0      =   lm/emax;
//     f0      =   x0*x0;
//     f0      =   f0*x0 - f0*f0/2;
//     el      =   max(root_finder_->FindRoot(
//                         x0, 1.0, 0.5,
//                         boost::bind(&Decay::Function, this, _1),
//                         boost::bind(&Decay::DifferentiatedFunction, this, _1),
//                         f0+(0.5-f0)*ernd) *emax
//                     , lm);
//     pl      =   sqrt(el*el - lm*lm);
//
//     if(store_neutrinos_)
//     {
//         int sign;
//         double cp, sp, ct, st, en, En1, En2;
//         cp  =   pl/(particle_->GetMass() - el);
//         sp  =   sqrt(1 - cp*cp);
//
//         if(arnd<0.5)
//         {
//             sign    =   1;
//             arnd    =   2*arnd;
//         }
//         else
//         {
//             sign    =   -1;
//             arnd    =   2*arnd-1;
//         }
//         ct  =   2*arnd - 1;
//         st  =   sqrt(1 - ct*ct);
//         en  =   (particle_->GetMass() - el)/2;
//         En1 =   -cp*ct - sign*sp*st;
//         En2 =   -cp*ct + sign*sp*st;
//         En1 =   en * ((particle_->GetEnergy()/particle_->GetMass()) + (particle_->GetMomentum()/particle_->GetMass()) * En1);
//         En2 =   en * ((particle_->GetEnergy()/particle_->GetMass()) + (particle_->GetMomentum()/particle_->GetMass()) * En2);
//         //o->output(1, out1, En1, 0);
//         //o->output(1, out2, En2, 0);
//     }
//
//     return el * (particle_->GetEnergy()/particle_->GetMass()) + pl * (particle_->GetMomentum()/particle_->GetMass()) * (2*arnd - 1);
// }


Decay::Decay():
    store_neutrinos_   ( false )
    ,multiplier_        ( 1 )
{
}

Decay::Decay(const Decay &decay)
    :store_neutrinos_   ( decay.store_neutrinos_ )
    ,multiplier_        ( decay.multiplier_ )
{
}

Decay& Decay::operator=(const Decay &decay)
{
    if (this != &decay)
    {
      Decay tmp(decay);
      swap(tmp);
    }
    return *this;
}

bool Decay::operator==(const Decay &decay) const
{
    if( store_neutrinos_    != decay.store_neutrinos_ ) return false;
    if( multiplier_         != decay.multiplier_ )      return false;

    //else
    return true;
}

bool Decay::operator!=(const Decay &decay) const
{
    return !(*this == decay);
}

void Decay::swap(Decay &decay)
{
    using std::swap;

    swap( store_neutrinos_  ,   decay.store_neutrinos_ );
    swap( multiplier_       ,   decay.multiplier_ );
}


// ------------------------------------------------------------------------- //
// Setter
// ------------------------------------------------------------------------- //

void Decay::SetStoreNeutrinos(bool store_neutrinos){
    store_neutrinos_    =   store_neutrinos;
}

void Decay::SetMultiplier(double multiplier){
    multiplier_ =   multiplier;
}

/******************************************************************************
*                               Decay Channels                                *
******************************************************************************/


bool DecayChannel::operator==(const DecayChannel& table) const
{
    return this->compare(table);
}

bool DecayChannel::operator!=(const DecayChannel& def) const
{
    return !(*this == def);
}

// DecayChannel::DecayChannel(const DecayChannel& mode)
// {
//     for (DecayProducts::iterator iter = decay_products.begin(); iter != decay_products.end(); ++iter)
//     {
//         delete *iter;
//     }
//
//     decay_products.clear();
//
//     for (DecayProducts::const_iterator iter = mode.decay_products.begin(); iter != mode.decay_products.end(); ++iter)
//     {
//         decay_products.push_back(new PROPOSALParticle(**iter));
//     }
// }

/******************************************************************************
*                                   Stable                                    *
******************************************************************************/

StableChannel::StableChannel()
    : DecayChannel()
{
}

StableChannel::~StableChannel()
{
}

StableChannel::StableChannel(const StableChannel& mode)
    : DecayChannel(mode)
{
}

bool StableChannel::compare(const DecayChannel& channel) const
{
    const StableChannel* stable = dynamic_cast<const StableChannel*>(&channel);

    if (!stable)
        return false;
    else
        return true;
}

DecayChannel::DecayProducts StableChannel::Decay(PROPOSALParticle*)
{
    // return empty vector;
    DecayProducts products;
    return products;
}

/******************************************************************************
*                                  Leponic                                   *
******************************************************************************/

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
    el      =   max(root_finder_.FindRoot(
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
    PROPOSALParticle* product_particle = new PROPOSALParticle(EMinusDef::Get());
    product_particle->SetEnergy(final_energy);
    product_particle->SetPosition(particle->GetPosition());
    product_particle->SetDirection(particle->GetDirection());
    product_particle->SetParticleId(particle->GetParticleId() + 1);
    product_particle->SetParentParticleId(particle->GetParentParticleId());
    product_particle->SetT(particle->GetT());
    product_particle->SetParentParticleEnergy(particle->GetEnergy());

    DecayProducts decay_products;
    decay_products.push_back(product_particle);

    return decay_products;
}

/******************************************************************************
*                             TwoBodyPhaseSpace                               *
******************************************************************************/

TwoBodyPhaseSpace::TwoBodyPhaseSpace(double m1, double m2)
    : first_daughter_mass_(m1)
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

DecayChannel::DecayProducts TwoBodyPhaseSpace::Decay(PROPOSALParticle* particle)
{
    double parent_mass = particle->GetMass();
    double el = (parent_mass * parent_mass + first_daughter_mass_ * first_daughter_mass_) / (2 * parent_mass);

    double arnd = RandomGenerator::Get().RandomDouble();

    double final_energy = el * (particle->GetEnergy() / particle->GetMass()) +
         sqrt(el * el - first_daughter_mass_ * first_daughter_mass_) * (particle->GetMomentum() / particle->GetMass()) * (2 * arnd - 1);

    // return el;
    // Create products
    PROPOSALParticle* product_particle = new PROPOSALParticle(EMinusDef::Get());
    product_particle->SetEnergy(final_energy);
    product_particle->SetPosition(particle->GetPosition());
    product_particle->SetDirection(particle->GetDirection());
    product_particle->SetParticleId(particle->GetParticleId() + 1);
    product_particle->SetParentParticleId(particle->GetParentParticleId());
    product_particle->SetT(particle->GetT());
    product_particle->SetParentParticleEnergy(particle->GetEnergy());

    DecayProducts decay_products;
    decay_products.push_back(product_particle);

    return decay_products;
}

