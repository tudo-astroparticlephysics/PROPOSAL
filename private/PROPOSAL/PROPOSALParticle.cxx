/*! \file   Particle.cxx
*   \brief  Source file for the Particle routines.
*
*   For more details see the class documentation.
*
*   \date   2013.03.14
*   \author Jan-Hendrik Köhne
*/


#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include "PROPOSAL/PROPOSALParticle.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Output.h"
#include <iomanip>


using namespace std;
using namespace PROPOSAL;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


PROPOSALParticle::PROPOSALParticle( )
    :propagated_distance_   ( 0 )
    ,x_                     ( 0 )
    ,y_                     ( 0 )
    ,z_                     ( 0 )
    ,t_                     ( 0 )
    ,theta_                 ( 0 )
    ,phi_                   ( 0 )
    ,costh_                 ( 1. )
    ,sinth_                 ( 0 )
    ,cosph_                 ( 1. )
    ,sinph_                 ( 0. )
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
    ,xi_                    ( 0 )
    ,yi_                    ( 0 )
    ,zi_                    ( 0 )
    ,ti_                    ( 0 )
    ,ei_                    ( 0 )
    ,xf_                    ( 0 )
    ,yf_                    ( 0 )
    ,zf_                    ( 0 )
    ,tf_                    ( 0 )
    ,ef_                    ( 0 )
    ,xc_                    ( 0 )
    ,yc_                    ( 0 )
    ,zc_                    ( 0 )
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
    ,x_                     ( particle.x_ )
    ,y_                     ( particle.y_ )
    ,z_                     ( particle.z_ )
    ,t_                     ( particle.t_ )
    ,theta_                 ( particle.theta_ )
    ,phi_                   ( particle.phi_ )
    ,costh_                 ( particle.costh_ )
    ,sinth_                 ( particle.sinth_ )
    ,cosph_                 ( particle.cosph_ )
    ,sinph_                 ( particle.sinph_ )
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
    ,xi_                    ( particle.xi_ )
    ,yi_                    ( particle.yi_ )
    ,zi_                    ( particle.zi_ )
    ,ti_                    ( particle.ti_ )
    ,ei_                    ( particle.ei_ )
    ,xf_                    ( particle.xf_ )
    ,yf_                    ( particle.yf_ )
    ,zf_                    ( particle.zf_ )
    ,tf_                    ( particle.tf_ )
    ,ef_                    ( particle.ef_ )
    ,xc_                    ( particle.xc_ )
    ,yc_                    ( particle.yc_ )
    ,zc_                    ( particle.zc_ )
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
                   double x,
                   double y,
                   double z,
                   double theta,
                   double phi,
                   double energy,
                   double t,
                   double prop_dist,
                   PROPOSALParticle *p)

    :propagated_distance_   ( prop_dist )
    ,x_                     ( x )
    ,y_                     ( y )
    ,z_                     ( z )
    ,t_                     ( t )
    ,theta_                 ( theta )
    ,phi_                   ( phi )
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
    ,xi_                    ( 0 )
    ,yi_                    ( 0 )
    ,zi_                    ( 0 )
    ,ti_                    ( 0 )
    ,ei_                    ( 0 )
    ,xf_                    ( 0 )
    ,yf_                    ( 0 )
    ,zf_                    ( 0 )
    ,tf_                    ( 0 )
    ,ef_                    ( 0 )
    ,xc_                    ( 0 )
    ,yc_                    ( 0 )
    ,zc_                    ( 0 )
    ,tc_                    ( 0 )
    ,ec_                    ( 0 )
    ,elost_                 ( 0 )
{
    InitParticle(type);
    SetEnergy(energy);
    Location(t, x, y, t, theta, phi);


    if(p!=NULL)
    {
        xi_      =   p->GetXi();
        yi_      =   p->GetYi();
        zi_      =   p->GetZi();
        ti_      =   p->GetTi();
        ei_      =   p->GetEi();
        xf_      =   p->GetXf();
        yf_      =   p->GetYf();
        zf_      =   p->GetZf();
        tf_      =   p->GetTf();
        ef_      =   p->GetEf();
        xc_      =   p->GetXc();
        yc_      =   p->GetYc();
        zc_      =   p->GetZc();
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
                   double x,
                   double y,
                   double z,
                   double theta,
                   double phi,
                   double energy,
                   double t,
                   double prop_dist,
                   double prim_energy)

    :propagated_distance_   ( prop_dist )
    ,x_                     ( x )
    ,y_                     ( y )
    ,z_                     ( z )
    ,t_                     ( t )
    ,theta_                 ( theta )
    ,phi_                   ( phi )
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
    ,xi_                    ( 0 )
    ,yi_                    ( 0 )
    ,zi_                    ( 0 )
    ,ti_                    ( 0 )
    ,ei_                    ( 0 )
    ,xf_                    ( 0 )
    ,yf_                    ( 0 )
    ,zf_                    ( 0 )
    ,tf_                    ( 0 )
    ,ef_                    ( 0 )
    ,xc_                    ( 0 )
    ,yc_                    ( 0 )
    ,zc_                    ( 0 )
    ,tc_                    ( 0 )
    ,ec_                    ( 0 )
    ,elost_                 ( 0 )

{
    InitParticle(type);
    SetEnergy(energy);
    Location(t, x, y, z, theta, phi);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


PROPOSALParticle::PROPOSALParticle(
                        ParticleType::Enum type,
                        double x,
                        double y,
                        double z,
                        double theta,
                        double phi,
                        double energy,
                        double t)

    :propagated_distance_   ( 0 )
    ,x_                     ( x )
    ,y_                     ( y )
    ,z_                     ( z )
    ,t_                     ( t )
    ,theta_                 ( theta )
    ,phi_                   ( phi )
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
    ,xi_                    ( 0 )
    ,yi_                    ( 0 )
    ,zi_                    ( 0 )
    ,ti_                    ( 0 )
    ,ei_                    ( 0 )
    ,xf_                    ( 0 )
    ,yf_                    ( 0 )
    ,zf_                    ( 0 )
    ,tf_                    ( 0 )
    ,ef_                    ( 0 )
    ,xc_                    ( 0 )
    ,yc_                    ( 0 )
    ,zc_                    ( 0 )
    ,tc_                    ( 0 )
    ,ec_                    ( 0 )
    ,elost_                 ( 0 )
{
    InitParticle(type);
    SetEnergy(energy);
    Location(t, x, y, z, theta, phi);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


PROPOSALParticle::PROPOSALParticle(ParticleType::Enum type)
    :propagated_distance_   ( 0 )
    ,x_                     ( 0 )
    ,y_                     ( 0 )
    ,z_                     ( 0 )
    ,t_                     ( 0 )
    ,theta_                 ( 0 )
    ,phi_                   ( 0 )
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
    ,xi_                    ( 0 )
    ,yi_                    ( 0 )
    ,zi_                    ( 0 )
    ,ti_                    ( 0 )
    ,ei_                    ( 0 )
    ,xf_                    ( 0 )
    ,yf_                    ( 0 )
    ,zf_                    ( 0 )
    ,tf_                    ( 0 )
    ,ef_                    ( 0 )
    ,xc_                    ( 0 )
    ,yc_                    ( 0 )
    ,zc_                    ( 0 )
    ,tc_                    ( 0 )
    ,ec_                    ( 0 )
    ,elost_                 ( 0 )
{
    InitParticle(type);
    SetEnergy(0);
    Location(0, 0, 0, 0, 0, 0);
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
    if(  x_                     != particle.x_)                     return false;
    if(  y_                     != particle.y_)                     return false;
    if(  z_                     != particle.z_)                     return false;
    if(  t_                     != particle.t_)                     return false;
    if(  theta_                 != particle.theta_)                 return false;
    if(  phi_                   != particle.phi_)                   return false;
    if(  costh_                 != particle.costh_)                 return false;
    if(  sinth_                 != particle.sinth_)                 return false;
    if(  cosph_                 != particle.cosph_)                 return false;
    if(  sinph_                 != particle.sinph_)                 return false;
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
    if(  xi_                    != particle.xi_)                    return false;
    if(  yi_                    != particle.yi_)                    return false;
    if(  zi_                    != particle.zi_)                    return false;
    if(  ti_                    != particle.ti_)                    return false;
    if(  ei_                    != particle.ei_)                    return false;
    if(  xf_                    != particle.xf_)                    return false;
    if(  yf_                    != particle.yf_)                    return false;
    if(  zf_                    != particle.zf_)                    return false;
    if(  tf_                    != particle.tf_)                    return false;
    if(  ef_                    != particle.ef_)                    return false;
    if(  xc_                    != particle.xc_)                    return false;
    if(  yc_                    != particle.yc_)                    return false;
    if(  zc_                    != particle.zc_)                    return false;
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
    os<<"\tx [cm]:\t\t\t\t\t"<<particle.x_<<endl;
    os<<"\ty [cm]:\t\t\t\t\t"<<particle.y_<<endl;
    os<<"\tz [cm]:\t\t\t\t\t"<<particle.z_<<endl;
    os<<"\ttheta [deg]:\t\t\t\t"<<particle.theta_<<endl;
    os<<"\tphi [deg]:\t\t\t\t"<<particle.phi_<<endl;
    os<<"\tage [s]:\t\t\t\t"<<particle.t_<<endl;
    os<<"\tparticle id:\t\t\t\t"<<particle.particle_id_<<endl;
    os<<"\tparent particle id_:\t\t\t"<<particle.parent_particle_id_<<scientific<<endl;
    os<<"\tenergy below paricle is lost [MeV]:\t"<<particle.low_<<fixed<<endl;
    os<<"\tpropagated distance [cm]:\t\t"<<particle.propagated_distance_<<endl;
    os<<"\tenergy lost in detector [MeV]:\t\t"<<particle.elost_<<endl;

    os<<"\n\tDetector entry point: x [cm] | y [cm] | z [cm] | time [s] | energy [MeV]"<<endl;
    os<<"\t\t"<<particle.xi_<<"\t"<<particle.yi_<<"\t"<<particle.zi_<<scientific<<"\t"<<particle.ti_<<fixed<<"\t"<<particle.ei_<<endl;
    os<<"\n\tDetector exit point: x [cm] | y [cm] | z [cm] | time [s] | energy [MeV]"<<endl;
    os<<"\t\t"<<particle.xf_<<"\t"<<particle.yf_<<"\t"<<particle.zf_<<scientific<<"\t"<<particle.tf_<<fixed<<"\t"<<particle.ef_<<endl;
    os<<"\n\tPoint of closest approach: x [cm] | y [cm] | z [cm] | time [s] | energy [MeV]"<<endl;
    os<<"\t\t"<<particle.xc_<<"\t"<<particle.yc_<<"\t"<<particle.zc_<<scientific<<"\t"<<particle.tc_<<fixed<<"\t"<<particle.ec_<<endl;
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
    swap( x_                     , particle.x_);
    swap( y_                     , particle.y_);
    swap( z_                     , particle.z_);
    swap( t_                     , particle.t_);
    swap( theta_                 , particle.theta_);
    swap( phi_                   , particle.phi_);
    swap( costh_                 , particle.costh_);
    swap( sinth_                 , particle.sinth_);
    swap( cosph_                 , particle.cosph_);
    swap( sinph_                 , particle.sinph_);
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
    swap( xi_                    , particle.xi_);
    swap( yi_                    , particle.yi_);
    swap( zi_                    , particle.zi_);
    swap( ti_                    , particle.ti_);
    swap( ei_                    , particle.ei_);
    swap( xf_                    , particle.xf_);
    swap( yf_                    , particle.yf_);
    swap( zf_                    , particle.zf_);
    swap( tf_                    , particle.tf_);
    swap( ef_                    , particle.ef_);
    swap( xc_                    , particle.xc_);
    swap( yc_                    , particle.yc_);
    swap( zc_                    , particle.zc_);
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
        default:
            log_warn("The particle type '%i' is set to Unknown!", type);
            type = ParticleType::unknown;
            name_ = "Unkown";
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
                        double x,
                        double y,
                        double z,
                        double theta,
                        double phi)
{
    //propagated_distance_ = 0; <--- seems to be wrong
    t_           =   time;
    x_           =   x;
    y_           =   y;
    z_           =   z;
    theta_       =   theta;
    phi_         =   phi;
    theta_       *=  (PI/180);
    phi_         *=  (PI/180);
    costh_       =   cos(theta_);
    sinth_       =   sin(theta_);
    cosph_       =   cos(phi_);
    sinph_       =   sin(phi_);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

void PROPOSALParticle::SetProperties(int parent_particle_id, int particle_id, double energy, double t,
                             double x, double y, double z, double theta, double phi,
                             double xi, double yi, double zi, double ti, double Ei,
                             double xf, double yf, double zf, double tf, double Ef,
                             double xc, double yc, double zc, double tc, double Ec)
{
    SetParentParticleId(parent_particle_id);
    SetParticleId(particle_id);
    SetEnergy(energy);
    SetT(t);
    SetX(x); SetY(y); SetZ(z); SetTheta(theta); SetPhi(phi);
    SetXi(xi); SetYi(yi); SetZi(zi); SetTi(ti); SetEi(Ei);
    SetXf(xf); SetYf(yf); SetZf(zf); SetTf(tf); SetEf(Ef);
    SetXc(xc); SetYc(yc); SetZc(zc); SetTc(tc); SetEc(Ec);
}

void PROPOSALParticle::SetEnergy(double e)
{
    energy_             =   e;

    if(energy_ < mass_) energy_ = mass_;

    square_momentum_    =   e*e-mass_*mass_;
    momentum_           =   sqrt(max(square_momentum_,0.0));
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void PROPOSALParticle::SetTheta(double theta)
{
    theta_  =   theta;
    theta_ *=   (PI/180);

    costh_  =   cos(theta_);
    sinth_  =   sin(theta_);

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void PROPOSALParticle::SetPhi(double phi)
{
    phi_   =    phi;
    phi_  *=    (PI/180);

    cosph_ =    cos(phi_);
    sinph_ =    sin(phi_);
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

void PROPOSALParticle::SetX(double x){
    x_ = x;
}

void PROPOSALParticle::SetY(double y){
    y_ = y;
}

void PROPOSALParticle::SetZ(double z){
    z_ = z;
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

void PROPOSALParticle::SetXi(double xi){
    xi_ = xi;
}

void PROPOSALParticle::SetYi(double yi){
    yi_ = yi;
}

void PROPOSALParticle::SetZi(double zi){
    zi_ = zi;
}

void PROPOSALParticle::SetTi(double ti){
    ti_ = ti;
}

void PROPOSALParticle::SetEi(double ei){
    ei_ = ei;
}

void PROPOSALParticle::SetXf(double xf){
    xf_ = xf;
}

void PROPOSALParticle::SetYf(double yf){
    yf_ = yf;
}

void PROPOSALParticle::SetZf(double zf){
    zf_ = zf;
}

void PROPOSALParticle::SetTf(double tf){
    tf_ = tf;
}

void PROPOSALParticle::SetEf(double ef){
    ef_ = ef;
}

void PROPOSALParticle::SetXc(double xc){
    xc_ = xc;
}

void PROPOSALParticle::SetYc(double yc){
    yc_ = yc;
}

void PROPOSALParticle::SetZc(double zc){
    zc_ = zc;
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
        if (particle_name.compare("mu-")) return ParticleType::MuMinus;
        else if (particle_name.compare("mu+")) return ParticleType::MuPlus;
        else return ParticleType::MuMinus;
    }
    else if (StartsWith(particle_name,"tau"))
    {
        if (particle_name.compare("tau-")) return ParticleType::TauMinus;
        else if (particle_name.compare("tau+")) return ParticleType::TauPlus;
        else return ParticleType::TauMinus;
    }
    else if (StartsWith(particle_name,"e"))
    {
        if (particle_name.compare("e-")) return ParticleType::EMinus;
        else if (particle_name.compare("e+")) return ParticleType::EPlus;
        else return ParticleType::EMinus;
    }
    else if (StartsWith(particle_name,"stau"))
    {
        if (particle_name.compare("stau-")) return ParticleType::STauMinus;
        else if (particle_name.compare("stau+")) return ParticleType::STauPlus;
        else return ParticleType::STauMinus;
    }
    else if (particle_name.compare("monopole")) return ParticleType::Monopole;
    else if (StartsWith(particle_name,"nu"))
    {
        if (particle_name.compare("nu_mu")) return ParticleType::NuMu;
        else if (particle_name.compare("nu_mu_bar")) return ParticleType::NuMuBar;
        else if (particle_name.compare("nu_tau")) return ParticleType::NuTau;
        else if (particle_name.compare("nu_tau_bar")) return ParticleType::NuTauBar;
        else if (particle_name.compare("nu_e")) return ParticleType::NuE;
        else if (particle_name.compare("nu_e_bar")) return ParticleType::NuEBar;
        else
        {
            log_fatal("the neutrino name '%s' is not correct it should be e.g. nu_mu_bar", particle_name.c_str());
        }
    }
    else if (particle_name.compare("DeltaE")) return ParticleType::DeltaE;
    else if (particle_name.compare("ContinuousEnergyLoss")) return ParticleType::ContinuousEnergyLoss;
    else if (particle_name.compare("Brems")) return ParticleType::Brems;
    else if (particle_name.compare("NuclInt")) return ParticleType::NuclInt;
    else if (particle_name.compare("Hadrons")) return ParticleType::Hadrons;
    else if (particle_name.compare("EPair")) return ParticleType::EPair;
    else
    {
        log_fatal("the particle name '%s' is not a PROPOSAL Particle", particle_name.c_str());
    }
}
