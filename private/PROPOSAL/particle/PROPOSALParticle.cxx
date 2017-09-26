/*! \file   Particle.cxx
*   \brief  Source file for the Particle routines.
*
*   For more details see the class documentation.
*
*   \date   2013.03.14
*   \author Jan-Hendrik Köhne
*/


// #include <cmath>
// #include <algorithm>
// #include <stdlib.h>
// #include <iomanip>

#include "PROPOSAL/particle/PROPOSALParticle.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Output.h"


using namespace std;
using namespace PROPOSAL;

/******************************************************************************
*                                  OStream                                    *
******************************************************************************/

namespace PROPOSAL
{

ostream& operator<<(ostream& os, PROPOSALParticle const& particle)
{
    os<<"---------------------------Particle( "<<&particle<<" )---------------------------"<<endl;
    os<<fixed<<setprecision(5);
    os<<"\tname:\t\t\t\t\t"<<particle.particle_def_.name<<endl;
    os<<"\tmass [MeV]:\t\t\t\t"<<particle.particle_def_.mass<<fixed<<endl;
    os<<"\tcharge:\t\t\t\t\t"<<particle.particle_def_.charge<<scientific<<endl;
    os<<"\tlifetime [s]:\t\t\t\t"<<particle.particle_def_.lifetime<<setprecision(8)<<endl;
    os<<"\tmomentum [MeV]:\t\t\t\t"<<particle.momentum_<<endl;
    os<<"\tenergy [MeV]: \t\t\t\t"<<particle.energy_<<fixed<<setprecision(5)<<endl;
    os<<"\tposition [cm]:\t\t\t\t\t"<<particle.position_<<endl;
    os<<"\tdirection :\t\t\t\t"<<particle.direction_<<endl;
    os<<"\tage [s]:\t\t\t\t"<<particle.time_<<endl;
    os<<"\tparticle id:\t\t\t\t"<<particle.particle_id_<<endl;
    os<<"\tparent particle id_:\t\t\t"<<particle.parent_particle_id_<<scientific<<endl;
    os<<"\tenergy below paricle is lost [MeV]:\t"<<particle.particle_def_.low<<fixed<<endl;
    os<<"\tpropagated distance [cm]:\t\t"<<particle.propagated_distance_<<endl;
    os<<"\tenergy lost in detector [MeV]:\t\t"<<particle.elost_<<endl;

    os<<"\n\tDetector entry point: x [cm] | y [cm] | z [cm] | time [s] | energy [MeV]"<<endl;
    os<<"\t\t"<<particle.entry_point_<<scientific<<"\t"<<particle.entry_time_<<fixed<<"\t"<<particle.entry_energy_<<endl;
    os<<"\n\tDetector exit point: x [cm] | y [cm] | z [cm] | time [s] | energy [MeV]"<<endl;
    os<<"\t\t"<<particle.exit_point_<<scientific<<"\t"<<particle.exit_time_<<fixed<<"\t"<<particle.exit_energy_<<endl;
    os<<"\n\tPoint of closest approach: x [cm] | y [cm] | z [cm] | time [s] | energy [MeV]"<<endl;
    os<<"\t\t"<<particle.closest_approach_point_<<scientific<<"\t"<<particle.closest_approach_time_<<fixed<<"\t"<<particle.closest_approach_energy_<<endl;
    os<<"--------------------------------------------------------------------------";
    return os;
}

}  // namespace PROPOSAL

/******************************************************************************
*                              Dynamic Particle                              *
******************************************************************************/

DynamicData::DynamicData(DynamicData::Type type)
    : type_id_(type)
    , position_(Vector3D())
    , direction_(Vector3D())
    , energy_(0)
    , parent_particle_energy_(0)
    , time_(0)
    , propagated_distance_(0)
{
}


DynamicData::DynamicData(const DynamicData& data)
    : type_id_(data.type_id_)
    , position_(data.position_)
    , direction_(data.direction_)
    , energy_(data.energy_)
    , parent_particle_energy_(data.parent_particle_energy_)
    , time_(data.time_)
    , propagated_distance_(data.propagated_distance_)
{
}

DynamicData::~DynamicData()
{
}


/******************************************************************************
*                              PROPOSALParticle                               *
******************************************************************************/

PROPOSALParticle::PROPOSALParticle()
    : DynamicData(DynamicData::Particle)
    , particle_def_(MuMinusDef::Get())
    , momentum_(0)
    , square_momentum_(0)
    , parent_particle_id_(0)
    , particle_id_(1)
    , entry_point_(Vector3D())
    , entry_time_(0)
    , entry_energy_(0)
    , exit_point_(Vector3D())
    , exit_time_(0)
    , exit_energy_(0)
    , closest_approach_point_(Vector3D())
    , closest_approach_time_(0)
    , closest_approach_energy_(0)
    , elost_(0)
{
    SetEnergy(energy_);
}

PROPOSALParticle::PROPOSALParticle(const ParticleDef& particleDef)
    : DynamicData(DynamicData::Particle)
    , particle_def_(particleDef)
    , momentum_(0)
    , square_momentum_(0)
    , parent_particle_id_(0)
    , particle_id_(1)
    , entry_point_(Vector3D())
    , entry_time_(0)
    , entry_energy_(0)
    , exit_point_(Vector3D())
    , exit_time_(0)
    , exit_energy_(0)
    , closest_approach_point_(Vector3D())
    , closest_approach_time_(0)
    , closest_approach_energy_(0)
    , elost_(0)
{
    SetEnergy(energy_);
}

// ------------------------------------------------------------------------- //
// Operators & swap
// ------------------------------------------------------------------------- //

// PROPOSALParticle& PROPOSALParticle::operator=(const PROPOSALParticle &particle){
//     if (this != &particle)
//     {
//       PROPOSALParticle tmp(particle);
//       swap(tmp);
//     }
//     return *this;
// }

// ------------------------------------------------------------------------- //
bool PROPOSALParticle::operator==(const PROPOSALParticle &particle) const
{
    if (propagated_distance_ != particle.propagated_distance_)
        return false;
    if (position_ != particle.position_)
        return false;
    if (direction_ != particle.direction_)
        return false;
    if (momentum_ != particle.momentum_)
        return false;
    if (square_momentum_ != particle.square_momentum_)
        return false;
    if (energy_ != particle.energy_)
        return false;
    if (parent_particle_id_ != particle.parent_particle_id_)
        return false;
    if (parent_particle_energy_ != particle.parent_particle_energy_)
        return false;
    if (particle_id_ != particle.particle_id_)
        return false;
    if (entry_point_ != particle.entry_point_)
        return false;
    if (entry_time_ != particle.entry_time_)
        return false;
    if (entry_energy_ != particle.entry_energy_)
        return false;
    if (exit_point_ != particle.exit_point_)
        return false;
    if (exit_time_ != particle.exit_time_)
        return false;
    if (exit_energy_ != particle.exit_energy_)
        return false;
    if (closest_approach_point_ != particle.closest_approach_point_)
        return false;
    if (closest_approach_time_ != particle.closest_approach_time_)
        return false;
    if (closest_approach_energy_ != particle.closest_approach_energy_)
        return false;
    if (elost_ != particle.elost_)
        return false;
    if (particle_def_ != particle.particle_def_)
    {
        return false;
    }

    // else
    return true;
}


// ------------------------------------------------------------------------- //
bool PROPOSALParticle::operator!=(const PROPOSALParticle &particle) const {
  return !(*this == particle);
}

// ------------------------------------------------------------------------- //
// void PROPOSALParticle::swap(PROPOSALParticle &particle)
// {
//     using std::swap;
//
//     swap( propagated_distance_   , particle.propagated_distance_);
//     position_.swap(particle.position_);
//     direction_.swap(particle.direction_);
//     swap( time_                  , particle.time_);
//     swap( momentum_              , particle.momentum_);
//     swap( square_momentum_       , particle.square_momentum_);
//     swap( energy_                , particle.energy_);
//     swap( parent_particle_id_    , particle.parent_particle_id_);
//     swap( parent_particle_energy_, particle.parent_particle_energy_);
//     swap( particle_id_           , particle.particle_id_);
//     entry_point_.swap(entry_point_);
//     swap( entry_time_            , particle.entry_time_);
//     swap( entry_energy_          , particle.entry_energy_);
//     exit_point_.swap(exit_point_);
//     swap( exit_time_             , particle.exit_time_);
//     swap( exit_energy_           , particle.exit_energy_);
//     closest_approach_point_.swap(closest_approach_point_);
//     swap( closest_approach_time_ , particle.closest_approach_time_);
//     swap( closest_approach_energy_ , particle.closest_approach_energy_);
//     swap( elost_                 , particle.elost_);
//     swap(particle_def_, particle.particle_def_);
// }


// ------------------------------------------------------------------------- //
// Setter
// ------------------------------------------------------------------------- //

void PROPOSALParticle::SetEnergy(double energy)
{
    energy_ = energy;

    if (energy_ < particle_def_.mass)
        energy_ = particle_def_.mass;

    square_momentum_ = energy * energy - particle_def_.mass * particle_def_.mass;
    momentum_ = sqrt(max(square_momentum_, 0.0));
}

// ------------------------------------------------------------------------- //
void PROPOSALParticle::SetMomentum(double momentum)
{
    momentum_ = momentum;
    square_momentum_ = momentum_ * momentum_;
    energy_ = sqrt(square_momentum_ + particle_def_.mass * particle_def_.mass);
}


// ParticleType::Enum PROPOSALParticle::GetTypeFromName(std::string particle_name)
// {
//     // returns the particle type of a particle name
//     // if there is just the name without charge (e.g. "mu" and not "mu+")
//     // then it is set to the default "mu-"
//
//     if (StartsWith(particle_name,"mu"))
//     {
//         if (particle_name.compare("mu-") == 0) return ParticleType::MuMinus;
//         else if (particle_name.compare("mu+") == 0) return ParticleType::MuPlus;
//         else return ParticleType::MuMinus;
//     }
//     else if (StartsWith(particle_name,"tau"))
//     {
//         if (particle_name.compare("tau-") == 0) return ParticleType::TauMinus;
//         else if (particle_name.compare("tau+") == 0) return ParticleType::TauPlus;
//         else return ParticleType::TauMinus;
//     }
//     else if (StartsWith(particle_name,"e"))
//     {
//         if (particle_name.compare("e-") == 0) return ParticleType::EMinus;
//         else if (particle_name.compare("e+") == 0) return ParticleType::EPlus;
//         else return ParticleType::EMinus;
//     }
//     else if (StartsWith(particle_name,"stau"))
//     {
//         if (particle_name.compare("stau-") == 0) return ParticleType::STauMinus;
//         else if (particle_name.compare("stau+") == 0) return ParticleType::STauPlus;
//         else return ParticleType::STauMinus;
//     }
//     else if (particle_name.compare("monopole") == 0) return ParticleType::Monopole;
//     else if (StartsWith(particle_name,"nu"))
//     {
//         if (particle_name.compare("nu_mu") == 0) return ParticleType::NuMu;
//         else if (particle_name.compare("nu_mu_bar") == 0) return ParticleType::NuMuBar;
//         else if (particle_name.compare("nu_tau") == 0) return ParticleType::NuTau;
//         else if (particle_name.compare("nu_tau_bar") == 0) return ParticleType::NuTauBar;
//         else if (particle_name.compare("nu_e") == 0) return ParticleType::NuE;
//         else if (particle_name.compare("nu_e_bar") == 0) return ParticleType::NuEBar;
//         else
//         {
//             log_fatal("the neutrino name '%s' is not correct it should be e.g. nu_mu_bar", particle_name.c_str());
//         }
//     }
//     else if (particle_name.compare("DeltaE") == 0) return ParticleType::DeltaE;
//     else if (particle_name.compare("ContinuousEnergyLoss") == 0) return ParticleType::ContinuousEnergyLoss;
//     else if (particle_name.compare("Brems") == 0) return ParticleType::Brems;
//     else if (particle_name.compare("NuclInt") == 0) return ParticleType::NuclInt;
//     else if (particle_name.compare("Hadrons") == 0) return ParticleType::Hadrons;
//     else if (particle_name.compare("EPair") == 0) return ParticleType::EPair;
//     else if (particle_name.compare("Unknown") == 0) return ParticleType::unknown;
//     else
//     {
//         log_fatal("the particle name '%s' is not a PROPOSAL Particle", particle_name.c_str());
//     }
// }
