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
    os<<"\tage [s]:\t\t\t\t"<<particle.t_<<endl;
    os<<"\tparticle id:\t\t\t\t"<<particle.particle_id_<<endl;
    os<<"\tparent particle id_:\t\t\t"<<particle.parent_particle_id_<<scientific<<endl;
    os<<"\tenergy below paricle is lost [MeV]:\t"<<particle.particle_def_.low<<fixed<<endl;
    os<<"\tpropagated distance [cm]:\t\t"<<particle.propagated_distance_<<endl;
    os<<"\tenergy lost in detector [MeV]:\t\t"<<particle.elost_<<endl;

    os<<"\n\tDetector entry point: x [cm] | y [cm] | z [cm] | time [s] | energy [MeV]"<<endl;
    os<<"\t\t"<<particle.entry_point_<<scientific<<"\t"<<particle.ti_<<fixed<<"\t"<<particle.ei_<<endl;
    os<<"\n\tDetector exit point: x [cm] | y [cm] | z [cm] | time [s] | energy [MeV]"<<endl;
    os<<"\t\t"<<particle.exit_point_<<scientific<<"\t"<<particle.tf_<<fixed<<"\t"<<particle.ef_<<endl;
    os<<"\n\tPoint of closest approach: x [cm] | y [cm] | z [cm] | time [s] | energy [MeV]"<<endl;
    os<<"\t\t"<<particle.closest_approach_point_<<scientific<<"\t"<<particle.tc_<<fixed<<"\t"<<particle.ec_<<endl;
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
    , t_(0)
    , propagated_distance_(0)
{
}


DynamicData::DynamicData(const DynamicData& data)
    : type_id_(data.type_id_)
    , position_(data.position_)
    , direction_(data.direction_)
    , energy_(data.energy_)
    , parent_particle_energy_(data.parent_particle_energy_)
    , t_(data.t_)
    , propagated_distance_(data.propagated_distance_)
{
}

DynamicData::~DynamicData()
{
}

// ------------------------------------------------------------------------- //
void DynamicData::SetPosition(const Vector3D& position)
{
    position_ = position;
}

// ------------------------------------------------------------------------- //
void DynamicData::SetDirection(const Vector3D& direction)
{
    direction_ = direction;
}

// ------------------------------------------------------------------------- //
void DynamicData::SetEnergy(double energy)
{
    energy_ = energy;
}

// ------------------------------------------------------------------------- //
void DynamicData::SetParentParticleEnergy(double parent_particle_energy){
    parent_particle_energy_ = parent_particle_energy;
}

// ------------------------------------------------------------------------- //
void DynamicData::SetT(double t){
    t_ = t;
}

// ------------------------------------------------------------------------- //
void DynamicData::SetPropagatedDistance(double prop_dist){
    propagated_distance_ = prop_dist;
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
    , ti_(0)
    , ei_(0)
    , exit_point_(Vector3D())
    , tf_(0)
    , ef_(0)
    , closest_approach_point_(Vector3D())
    , tc_(0)
    , ec_(0)
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
    , ti_(0)
    , ei_(0)
    , exit_point_(Vector3D())
    , tf_(0)
    , ef_(0)
    , closest_approach_point_(Vector3D())
    , tc_(0)
    , ec_(0)
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
    if (t_ != particle.t_)
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
    if (ti_ != particle.ti_)
        return false;
    if (ei_ != particle.ei_)
        return false;
    if (exit_point_ != particle.exit_point_)
        return false;
    if (tf_ != particle.tf_)
        return false;
    if (ef_ != particle.ef_)
        return false;
    if (closest_approach_point_ != particle.closest_approach_point_)
        return false;
    if (tc_ != particle.tc_)
        return false;
    if (ec_ != particle.ec_)
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
//     swap( t_                     , particle.t_);
//     direction_.swap(particle.direction_);
//     swap( momentum_              , particle.momentum_);
//     swap( square_momentum_       , particle.square_momentum_);
//     swap( energy_                , particle.energy_);
//     swap( parent_particle_id_    , particle.parent_particle_id_);
//     swap( parent_particle_energy_, particle.parent_particle_energy_);
//     swap( particle_id_           , particle.particle_id_);
//     entry_point_.swap(entry_point_);
//     swap( ti_                    , particle.ti_);
//     swap( ei_                    , particle.ei_);
//     exit_point_.swap(exit_point_);
//     swap( tf_                    , particle.tf_);
//     swap( ef_                    , particle.ef_);
//     closest_approach_point_.swap(closest_approach_point_);
//     swap( tc_                    , particle.tc_);
//     swap( ec_                    , particle.ec_);
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

// ------------------------------------------------------------------------- //
void PROPOSALParticle::SetParentParticleId(int parent_particle_id){
    parent_particle_id_ = parent_particle_id;
}

// ------------------------------------------------------------------------- //
void PROPOSALParticle::SetParticleId(int particle_id){
    particle_id_ = particle_id;
}

// ------------------------------------------------------------------------- //
void PROPOSALParticle::SetEntryPoint(Vector3D& entry_point)
{
    entry_point_ = entry_point;
}

// ------------------------------------------------------------------------- //
void PROPOSALParticle::SetTi(double ti){
    ti_ = ti;
}

// ------------------------------------------------------------------------- //
void PROPOSALParticle::SetEi(double ei){
    ei_ = ei;
}

// ------------------------------------------------------------------------- //
void PROPOSALParticle::SetExitPoint(Vector3D& exit_point)
{
    exit_point_ = exit_point;
}

// ------------------------------------------------------------------------- //
void PROPOSALParticle::SetTf(double tf){
    tf_ = tf;
}

// ------------------------------------------------------------------------- //
void PROPOSALParticle::SetEf(double ef){
    ef_ = ef;
}

// ------------------------------------------------------------------------- //
void PROPOSALParticle::SetClosestApproachPoint(Vector3D& closest_approach_point)
{
    closest_approach_point_ = closest_approach_point;
}

// ------------------------------------------------------------------------- //
void PROPOSALParticle::SetTc(double tc){
    tc_ = tc;
}

// ------------------------------------------------------------------------- //
void PROPOSALParticle::SetEc(double ec){
    ec_ = ec;
}

// ------------------------------------------------------------------------- //
void PROPOSALParticle::SetElost(double elost){
    elost_ = elost;
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
