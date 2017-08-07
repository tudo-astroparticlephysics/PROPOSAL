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

#include "PROPOSAL/PROPOSALParticle.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Output.h"


using namespace std;
using namespace PROPOSAL;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


PROPOSALParticle::PROPOSALParticle( )
    :propagated_distance_   ( 0 )
    ,position_              ( Vector3D() )
    ,t_                     ( 0 )
    ,direction_             ( Vector3D() )
    ,momentum_              ( 0 )
    ,square_momentum_       ( 0 )
    ,energy_                ( 0 )
    ,mass_                  ( MMU )
    ,lifetime_              ( LMU )
    ,charge_                ( -1 )
    ,name_                  ( "mu-" )
    ,low_                   ( mass_ )
    ,type_                  ( ParticleType::MuMinus )
    ,parent_particle_id_    ( 0 )
    ,parent_particle_energy_( 0 )
    ,particle_id_           ( 1 )
    ,entry_point_           ( Vector3D() )
    ,ti_                    ( 0 )
    ,ei_                    ( 0 )
    ,exit_point_            ( Vector3D() )
    ,tf_                    ( 0 )
    ,ef_                    ( 0 )
    ,closest_approach_point_ ( Vector3D() )
    ,tc_                    ( 0 )
    ,ec_                    ( 0 )
    ,elost_                 ( 0 )
{
    SetEnergy( energy_ );
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


PROPOSALParticle::PROPOSALParticle(const PROPOSALParticle& particle)
    :propagated_distance_   ( particle.propagated_distance_ )
    ,position_              ( particle.position_ )
    ,t_                     ( particle.t_ )
    ,direction_             ( particle.direction_ )
    ,momentum_              ( particle.momentum_ )
    ,square_momentum_       ( particle.square_momentum_ )
    ,energy_                ( particle.energy_ )
    ,mass_                  ( particle.mass_ )
    ,lifetime_              ( particle.lifetime_ )
    ,charge_                ( particle.charge_ )
    ,name_                  ( particle.name_ )
    ,low_                   ( particle.low_ )
    ,type_                  ( particle.type_ )
    ,parent_particle_id_    ( particle.parent_particle_id_ )
    ,parent_particle_energy_( particle.parent_particle_energy_  )
    ,particle_id_           ( particle.particle_id_ )
    ,entry_point_           ( particle.entry_point_ )
    ,ti_                    ( particle.ti_ )
    ,ei_                    ( particle.ei_ )
    ,exit_point_            ( particle.exit_point_ )
    ,tf_                    ( particle.tf_ )
    ,ef_                    ( particle.ef_ )
    ,closest_approach_point_ ( particle.closest_approach_point_ )
    ,tc_                    ( particle.tc_ )
    ,ec_                    ( particle.ec_ )
    ,elost_                 ( particle.elost_ )
{

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


PROPOSALParticle::PROPOSALParticle(int parent_particle_id,
                   int particle_id,
                   ParticleType::Enum type,
                   Vector3D position,
                   Vector3D direction,
                   double energy,
                   double t,
                   double prop_dist,
                   PROPOSALParticle *p)

    :propagated_distance_   ( prop_dist )
    ,position_              ( position )
    ,t_                     ( t )
    ,direction_             (direction)
    ,momentum_              ( 0 )
    ,square_momentum_       ( 0 )
    ,energy_                ( 0 )
    ,mass_                  ( MMU )
    ,lifetime_              ( LMU )
    ,charge_                ( -1 )
    ,name_                  ( "mu-" )
    ,low_                   ( mass_ )
    ,type_                  ( ParticleType::MuMinus )
    ,parent_particle_id_    ( parent_particle_id )
    ,parent_particle_energy_( 0 )
    ,particle_id_           ( particle_id )
    ,entry_point_           ( Vector3D() )
    ,ti_                    ( 0 )
    ,ei_                    ( 0 )
    ,exit_point_            ( Vector3D() )
    ,tf_                    ( 0 )
    ,ef_                    ( 0 )
    ,closest_approach_point_ ( Vector3D() )
    ,tc_                    ( 0 )
    ,ec_                    ( 0 )
    ,elost_                 ( 0 )
{
    InitParticle(type);
    SetEnergy(energy);
    Location(t, position, direction);


    if(p!=NULL)
    {
        entry_point_ = p->GetEntryPoint();
        ti_      =   p->GetTi();
        ei_      =   p->GetEi();
        exit_point_ = p->GetExitPoint();
        tf_      =   p->GetTf();
        ef_      =   p->GetEf();
        closest_approach_point_ = p->GetClosestApproachPoint();
        tc_      =   p->GetTc();
        ec_      =   p->GetEc();
        elost_   =   p->GetElost();
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


PROPOSALParticle::PROPOSALParticle(int parent_particle_id,
                   int particle_id,
                   ParticleType::Enum type,
                   Vector3D position,
                   Vector3D direction,
                   double energy,
                   double t,
                   double prop_dist,
                   double prim_energy)

    :propagated_distance_   ( prop_dist )
    ,position_              ( position )
    ,t_                     ( t )
    ,direction_             (direction)
    ,momentum_              ( 0 )
    ,square_momentum_       ( 0 )
    ,energy_                ( energy )
    ,mass_                  ( MMU )
    ,lifetime_              ( LMU )
    ,charge_                ( -1 )
    ,name_                  ( "mu-" )
    ,low_                   ( mass_ )
    ,type_                  ( ParticleType::MuMinus )
    ,parent_particle_id_    ( parent_particle_id )
    ,parent_particle_energy_( prim_energy )
    ,particle_id_           ( particle_id )
    ,entry_point_           ( Vector3D() )
    ,ti_                    ( 0 )
    ,ei_                    ( 0 )
    ,exit_point_            ( Vector3D() )
    ,tf_                    ( 0 )
    ,ef_                    ( 0 )
    ,closest_approach_point_ ( Vector3D() )
    ,tc_                    ( 0 )
    ,ec_                    ( 0 )
    ,elost_                 ( 0 )

{
    InitParticle(type);
    SetEnergy(energy);
    Location(t, position, direction);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


PROPOSALParticle::PROPOSALParticle(
                        ParticleType::Enum type,
                        Vector3D position,
                        Vector3D direction,
                        double energy,
                        double t)

    :propagated_distance_   ( 0 )
    ,position_              ( position )
    ,t_                     ( t )
    ,direction_             (direction)
    ,momentum_              ( 0 )
    ,square_momentum_       ( 0 )
    ,energy_                ( energy )
    ,mass_                  ( MMU )
    ,lifetime_              ( LMU )
    ,charge_                ( -1 )
    ,name_                  ( "mu-" )
    ,low_                   ( mass_ )
    ,type_                  ( ParticleType::MuMinus )
    ,parent_particle_id_    ( 0 )
    ,parent_particle_energy_( 0 )
    ,particle_id_           ( 1 )
    ,entry_point_           ( Vector3D() )
    ,ti_                    ( 0 )
    ,ei_                    ( 0 )
    ,exit_point_            ( Vector3D() )
    ,tf_                    ( 0 )
    ,ef_                    ( 0 )
    ,closest_approach_point_ ( Vector3D() )
    ,tc_                    ( 0 )
    ,ec_                    ( 0 )
    ,elost_                 ( 0 )
{
    InitParticle(type);
    SetEnergy(energy);
    Location(t, position, direction);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


PROPOSALParticle::PROPOSALParticle(ParticleType::Enum type)
    :propagated_distance_   ( 0 )
    ,position_              ( Vector3D() )
    ,t_                     ( 0 )
    ,direction_             ( Vector3D() )
    ,momentum_              ( 0 )
    ,square_momentum_       ( 0 )
    ,energy_                ( 0 )
    ,mass_                  ( MMU )
    ,lifetime_              ( LMU )
    ,charge_                ( -1 )
    ,name_                  ( "mu-" )
    ,low_                   ( mass_ )
    ,type_                  ( ParticleType::MuMinus )
    ,parent_particle_id_    ( 0 )
    ,parent_particle_energy_( 0 )
    ,particle_id_           ( 1 )
    ,entry_point_           ( Vector3D() )
    ,ti_                    ( 0 )
    ,ei_                    ( 0 )
    ,exit_point_            ( Vector3D() )
    ,tf_                    ( 0 )
    ,ef_                    ( 0 )
    ,closest_approach_point_ ( Vector3D() )
    ,tc_                    ( 0 )
    ,ec_                    ( 0 )
    ,elost_                 ( 0 )
{
    InitParticle(type);
    SetEnergy(0);
    Location(0, position_, direction_);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


PROPOSALParticle& PROPOSALParticle::operator=(const PROPOSALParticle &particle){
    if (this != &particle)
    {
      PROPOSALParticle tmp(particle);
      swap(tmp);
    }
    return *this;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool PROPOSALParticle::operator==(const PROPOSALParticle &particle) const
{
    if(  propagated_distance_   != particle.propagated_distance_)   return false;
    if(  position_              != particle.position_)              return false;
    if(  t_                     != particle.t_)                     return false;
    if(  direction_             != particle.direction_)             return false;
    if(  momentum_              != particle.momentum_)              return false;
    if(  square_momentum_       != particle.square_momentum_)       return false;
    if(  energy_                != particle.energy_)                return false;
    if(  mass_                  != particle.mass_)                  return false;
    if(  lifetime_              != particle.lifetime_)              return false;
    if(  charge_                != particle.charge_)                return false;
    if(  low_                   != particle.low_)                   return false;
    if(  type_                  != particle.type_)                  return false;
    if(  parent_particle_id_    != particle.parent_particle_id_)    return false;
    if(  parent_particle_energy_!= particle.parent_particle_energy_)return false;
    if(  particle_id_           != particle.particle_id_)           return false;
    if(  entry_point_           != particle.entry_point_)           return false;
    if(  ti_                    != particle.ti_)                    return false;
    if(  ei_                    != particle.ei_)                    return false;
    if(  exit_point_            != particle.exit_point_)            return false;
    if(  tf_                    != particle.tf_)                    return false;
    if(  ef_                    != particle.ef_)                    return false;
    if(  closest_approach_point_!= particle.closest_approach_point_)return false;
    if(  tc_                    != particle.tc_)                    return false;
    if(  ec_                    != particle.ec_)                    return false;
    if(  elost_                 != particle.elost_)                 return false;

    if(  name_.compare(particle.name_) != 0)        return false;

    //else
    return true;

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool PROPOSALParticle::operator!=(const PROPOSALParticle &particle) const {
  return !(*this == particle);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

namespace PROPOSAL
{

ostream& operator<<(ostream& os, PROPOSALParticle const& particle)
{
    os<<"---------------------------Particle( "<<&particle<<" )---------------------------"<<endl;
    os<<fixed<<setprecision(5);
    os<<"\tname:\t\t\t\t\t"<<particle.name_<<endl;
    os<<"\ttype:\t\t\t\t\t"<<particle.type_<<scientific<<endl;
    os<<"\tmass [MeV]:\t\t\t\t"<<particle.mass_<<fixed<<endl;
    os<<"\tcharge:\t\t\t\t\t"<<particle.charge_<<scientific<<endl;
    os<<"\tlifetime [s]:\t\t\t\t"<<particle.lifetime_<<setprecision(8)<<endl;
    os<<"\tmomentum [MeV]:\t\t\t\t"<<particle.momentum_<<endl;
    os<<"\tenergy [MeV]: \t\t\t\t"<<particle.energy_<<fixed<<setprecision(5)<<endl;
    os<<"\tposition [cm]:\t\t\t\t\t"<<particle.position_<<endl;
    os<<"\tdirection :\t\t\t\t"<<particle.direction_<<endl;
    os<<"\tage [s]:\t\t\t\t"<<particle.t_<<endl;
    os<<"\tparticle id:\t\t\t\t"<<particle.particle_id_<<endl;
    os<<"\tparent particle id_:\t\t\t"<<particle.parent_particle_id_<<scientific<<endl;
    os<<"\tenergy below paricle is lost [MeV]:\t"<<particle.low_<<fixed<<endl;
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

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void PROPOSALParticle::swap(PROPOSALParticle &particle)
{
    using std::swap;

    swap( propagated_distance_   , particle.propagated_distance_);
    position_.swap(particle.position_);
    swap( t_                     , particle.t_);
    direction_.swap(particle.direction_);
    swap( momentum_              , particle.momentum_);
    swap( square_momentum_       , particle.square_momentum_);
    swap( energy_                , particle.energy_);
    swap( mass_                  , particle.mass_);
    swap( lifetime_              , particle.lifetime_);
    swap( charge_                , particle.charge_);
    swap( low_                   , particle.low_);
    swap( type_                  , particle.type_);
    swap( parent_particle_id_    , particle.parent_particle_id_);
    swap( parent_particle_energy_, particle.parent_particle_energy_);
    swap( particle_id_           , particle.particle_id_);
    entry_point_.swap(entry_point_);
    swap( ti_                    , particle.ti_);
    swap( ei_                    , particle.ei_);
    exit_point_.swap(exit_point_);
    swap( tf_                    , particle.tf_);
    swap( ef_                    , particle.ef_);
    closest_approach_point_.swap(closest_approach_point_);
    swap( tc_                    , particle.tc_);
    swap( ec_                    , particle.ec_);
    swap( elost_                 , particle.elost_);

    name_.swap(particle.name_);


}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//------------------------private member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void PROPOSALParticle::InitParticle(ParticleType::Enum type){

    switch (type)
    {
        case ParticleType::MuMinus:
            name_ = "mu-";
            mass_       =   MMU;
            lifetime_   =   LMU;
            charge_     =   -1;
            break;
        case ParticleType::MuPlus:
            name_ = "mu+";
            mass_       =   MMU;
            lifetime_   =   LMU;
            charge_     =   1;
            break;
        case ParticleType::EMinus:
            name_ = "e-";
            mass_       =   ME;
            lifetime_   =   -1;
            charge_     =   -1;
            break;
        case ParticleType::EPlus:
            name_ = "e+";
            mass_       =   ME;
            lifetime_   =   -1;
            charge_     =   1;
            break;
        case ParticleType::TauMinus:
            name_ = "tau-";
            mass_       =   MTAU;
            lifetime_   =   LTAU;
            charge_     =   -1;
            break;
        case ParticleType::TauPlus:
            name_ = "tau+";
            mass_       =   MTAU;
            lifetime_   =   LTAU;
            charge_     =   1;
            break;
        case ParticleType::STauMinus:
            name_ = "stau-";
            mass_       =   MSTAU;
            lifetime_   =   LSTAU;
            charge_     =   -1;
            break;
        case ParticleType::STauPlus:
            name_ = "stau+";
            mass_       =   MSTAU;
            lifetime_   =   LSTAU;
            charge_     =   1;
            break;
        case ParticleType::NuE:
            name_ = "nu_e";
            mass_       =   0;
            lifetime_   =   -1;
            charge_     =   0;
            break;
        case ParticleType::NuEBar:
            name_ = "~nu_e";
            mass_       =   0;
            lifetime_   =   -1;
            charge_     =   0;
            break;
        case ParticleType::NuMu:
            name_ = "nu_mu";
            mass_       =   0;
            lifetime_   =   -1;
            charge_     =   0;
            break;
        case ParticleType::NuMuBar:
            name_ = "~nu_mu";
            mass_       =   0;
            lifetime_   =   -1;
            charge_     =   0;
            break;
        case ParticleType::NuTau:
            name_ = "nu_tau";
            mass_       =   0;
            lifetime_   =   -1;
            charge_     =   0;
            break;
        case ParticleType::NuTauBar:
            name_ = "~nu_tau";
            mass_       =   0;
            lifetime_   =   -1;
            charge_     =   0;
            break;
        case ParticleType::DeltaE:
            name_ = "DeltaE";
            mass_       =   0;
            lifetime_   =   0;
            break;
        case ParticleType::Brems:
            name_ = "Brems";
            mass_       =   0;
            lifetime_   =   0;
            break;
        case ParticleType::NuclInt:
            name_ = "NuclInt";
            mass_       =   0;
            lifetime_   =   0;
            break;
        case ParticleType::Hadrons:
            name_ = "Hadrons";
            mass_       =   0;
            lifetime_   =   0;
            break;
        case ParticleType::ContinuousEnergyLoss:
            name_ = "ContinuousEnergyLoss";
            mass_       =   0;
            lifetime_   =   0;
            break;
        case ParticleType::EPair:
            name_ = "EPair";
            mass_       =   0;
            lifetime_   =   0;
            break;
        case ParticleType::MuPair:
            name_ = "MuPair";
            mass_       =   0;
            lifetime_   =   0;
            break;
        case ParticleType::Monopole:
            name_ = "monopole";
            mass_       =   MMON;
            lifetime_   =   -1;
            charge_     =   CMON;
            break;
        case ParticleType::Gamma:
            name_ = "Gamma";
            mass_       =   0;
            lifetime_   =   -1;
            charge_     =   0;
            break;
        case ParticleType::Pi0:
            name_ = "Pi0";
            mass_       =   MPI0;
            lifetime_   =   LPI0;
            charge_     =   0;
            break;
        case ParticleType::PiPlus:
            name_ = "PiPlus";
            mass_       =   MPI;
            lifetime_   =   LPI;
            charge_     =   1;
            break;
        case ParticleType::PiMinus:
            name_ = "PiMinus";
            mass_       =   MPI;
            lifetime_   =   LPI;
            charge_     =   -1;
            break;
        case ParticleType::KPlus:
            name_ = "KPlus";
            mass_       =   MKAON;
            lifetime_   =   LKAON;
            charge_     =   1;
            break;
        case ParticleType::KMinus:
            name_ = "KMinus";
            mass_       =   MKAON;
            lifetime_   =   LKAON;
            charge_     =   -1;
            break;
        case ParticleType::PPlus:
            name_ = "PPlus";
            mass_       =   MP;
            lifetime_   =   -1;
            charge_     =   1;
            break;
        case ParticleType::PMinus:
            name_ = "PMinus";
            mass_       =   MP;
            lifetime_   =   -1;
            charge_     =   -1;
            break;
        case ParticleType::StableMassiveParticle:
            name_ = "SMP";
            mass_       =   MSMP;
            lifetime_   =   LSMP;
            charge_     =   -1;
            break;
        default:
            log_warn("The particle type '%i' is set to Unknown!", type);
            type = ParticleType::unknown;
            name_ = "Unknown";
            mass_ = 0;
            lifetime_ = 0;
            break;
    }

    type_ = type;
    low_ = mass_;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void PROPOSALParticle::Location(double time,
                        Vector3D& position,
                        Vector3D& direction)
{
    //propagated_distance_ = 0; <--- seems to be wrong
    t_           =   time;
    position_    =   position;
    direction_   =   direction;

    // costh_       =   cos(direction.GetTheta());
    // sinth_       =   sin(direction.GetTheta());
    // cosph_       =   cos(direction.GetPhi());
    // sinph_       =   sin(direction.GetPhi());
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

// void PROPOSALParticle::SetProperties(int parent_particle_id, int particle_id, double energy, double t,
//                                     Vector3D& position, Vector3D& direction,
//                                     Vector3D& entry_point, double ti, double Ei,
//                                     Vector3D& exit_point, double tf, double Ef,
//                                     Vector3D& closest_approach_point, double tc, double Ec)
// {
//     SetParentParticleId(parent_particle_id);
//     SetParticleId(particle_id);
//     SetEnergy(energy);
//     SetT(t);
//     SetPosition(position);
//     SetDirection(direction);
//     SetEntryPoint(entry_point); SetTi(ti); SetEi(Ei);
//     SetExitPoint(exit_point); SetTf(tf); SetEf(Ef);
//     SetClosestApproachPoint(closest_approach_point); SetTc(tc); SetEc(Ec);
// }

void PROPOSALParticle::SetEnergy(double energy)
{
    energy_ = energy;

    if(energy_ < mass_) energy_ = mass_;

    square_momentum_ = energy*energy-mass_*mass_;
    momentum_ = sqrt(max(square_momentum_,0.0));
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void PROPOSALParticle::SetDirection(Vector3D& direction)
{
    direction_ = direction;

    // costh_ = cos(direction.GetTheta());
    // sinth_ = sin(direction.GetTheta());

    // cosph_ = cos(direction.GetPhi());
    // sinph_ = sin(direction.GetPhi());
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void PROPOSALParticle::SetMomentum(double momentum)
{
    momentum_           =   momentum;
    square_momentum_    =   momentum_*momentum_;
    energy_             =   sqrt(square_momentum_ + mass_*mass_);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void PROPOSALParticle::SetPropagatedDistance(double prop_dist){
    propagated_distance_ = prop_dist;
}

void PROPOSALParticle::SetPosition(Vector3D& position)
{
    position_ = position;
}

void PROPOSALParticle::SetT(double t){
    t_ = t;
}

void PROPOSALParticle::SetMass(double mass){
    mass_ = mass;

    if (low_ < mass)
    {
        low_ = mass;
    }
    else
    {
        log_warn("Your lower interpolation bound is higher than your mass: %f > %f", low_, mass_);
    }
}

void PROPOSALParticle::SetLifetime(double lifetime){
    lifetime_ = lifetime;
}

void PROPOSALParticle::SetCharge(double charge){
    charge_ = charge;
}

void PROPOSALParticle::SetName(std::string name){
    name_ = name;
}

void PROPOSALParticle::SetLow(double low){
    low_ = low;
}

void PROPOSALParticle::SetType(ParticleType::Enum type){
    type_ = type;
}

void PROPOSALParticle::SetParentParticleId(int parent_particle_id){
    parent_particle_id_ = parent_particle_id;
}

void PROPOSALParticle::SetParentParticleEnergy(double parent_particle_energy){
    parent_particle_energy_ = parent_particle_energy;
}

void PROPOSALParticle::SetParticleId(int particle_id){
    particle_id_ = particle_id;
}

void PROPOSALParticle::SetEntryPoint(Vector3D& entry_point)
{
    entry_point_ = entry_point;
}

void PROPOSALParticle::SetTi(double ti){
    ti_ = ti;
}

void PROPOSALParticle::SetEi(double ei){
    ei_ = ei;
}

void PROPOSALParticle::SetExitPoint(Vector3D& exit_point)
{
    exit_point_ = exit_point;
}

void PROPOSALParticle::SetTf(double tf){
    tf_ = tf;
}

void PROPOSALParticle::SetEf(double ef){
    ef_ = ef;
}

void PROPOSALParticle::SetClosestApproachPoint(Vector3D& closest_approach_point)
{
    closest_approach_point_ = closest_approach_point;
}

void PROPOSALParticle::SetTc(double tc){
    tc_ = tc;
}

void PROPOSALParticle::SetEc(double ec){
    ec_ = ec;
}

void PROPOSALParticle::SetElost(double elost){
    elost_ = elost;
}

std::string PROPOSALParticle::GetName(ParticleType::Enum pt) {

    PROPOSALParticle p(pt);
    return p.GetName();
}

ParticleType::Enum PROPOSALParticle::GetTypeFromName(std::string particle_name)
{
    // returns the particle type of a particle name
    // if there is just the name without charge (e.g. "mu" and not "mu+")
    // then it is set to the default "mu-"

    if (StartsWith(particle_name,"mu"))
    {
        if (particle_name.compare("mu-") == 0) return ParticleType::MuMinus;
        else if (particle_name.compare("mu+") == 0) return ParticleType::MuPlus;
        else return ParticleType::MuMinus;
    }
    else if (StartsWith(particle_name,"tau"))
    {
        if (particle_name.compare("tau-") == 0) return ParticleType::TauMinus;
        else if (particle_name.compare("tau+") == 0) return ParticleType::TauPlus;
        else return ParticleType::TauMinus;
    }
    else if (StartsWith(particle_name,"e"))
    {
        if (particle_name.compare("e-") == 0) return ParticleType::EMinus;
        else if (particle_name.compare("e+") == 0) return ParticleType::EPlus;
        else return ParticleType::EMinus;
    }
    else if (StartsWith(particle_name,"stau"))
    {
        if (particle_name.compare("stau-") == 0) return ParticleType::STauMinus;
        else if (particle_name.compare("stau+") == 0) return ParticleType::STauPlus;
        else return ParticleType::STauMinus;
    }
    else if (particle_name.compare("monopole") == 0) return ParticleType::Monopole;
    else if (StartsWith(particle_name,"nu"))
    {
        if (particle_name.compare("nu_mu") == 0) return ParticleType::NuMu;
        else if (particle_name.compare("nu_mu_bar") == 0) return ParticleType::NuMuBar;
        else if (particle_name.compare("nu_tau") == 0) return ParticleType::NuTau;
        else if (particle_name.compare("nu_tau_bar") == 0) return ParticleType::NuTauBar;
        else if (particle_name.compare("nu_e") == 0) return ParticleType::NuE;
        else if (particle_name.compare("nu_e_bar") == 0) return ParticleType::NuEBar;
        else
        {
            log_fatal("the neutrino name '%s' is not correct it should be e.g. nu_mu_bar", particle_name.c_str());
        }
    }
    else if (particle_name.compare("DeltaE") == 0) return ParticleType::DeltaE;
    else if (particle_name.compare("ContinuousEnergyLoss") == 0) return ParticleType::ContinuousEnergyLoss;
    else if (particle_name.compare("Brems") == 0) return ParticleType::Brems;
    else if (particle_name.compare("NuclInt") == 0) return ParticleType::NuclInt;
    else if (particle_name.compare("Hadrons") == 0) return ParticleType::Hadrons;
    else if (particle_name.compare("EPair") == 0) return ParticleType::EPair;
    else if (particle_name.compare("Unknown") == 0) return ParticleType::unknown;
    else
    {
        log_fatal("the particle name '%s' is not a PROPOSAL Particle", particle_name.c_str());
    }
}
