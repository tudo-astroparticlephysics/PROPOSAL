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
#include "PROPOSAL/Particle.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Constants.h"
#include <iomanip>


using namespace std;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Particle::Particle( )
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
    ,charge_                ( 1 )
    ,name_                  ( "mu" )
    ,low_                   ( mass_ )
    ,type_                  ( 1 )
    ,parent_particle_id_    ( 0 )
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


Particle::Particle(const Particle& particle)
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


Particle::Particle(int parent_particle_id,
                   int particle_id,
                   string name,
                   double x,
                   double y,
                   double z,
                   double theta,
                   double phi,
                   double energy,
                   double t,
                   double prop_dist,
                   Particle *p)

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
    ,charge_                ( 1 )
    ,name_                  ( "mu" )
    ,low_                   ( mass_ )
    ,type_                  ( 1 )
    ,parent_particle_id_    ( parent_particle_id )
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
    if(StartsWith(name,"tau"))
    {
        name_        =   "tau";
        type_        =   2;
        mass_        =   MTAU;
        lifetime_    =   LTAU;
    }
    else if(StartsWith(name,"e"))
    {
        name_        =   "e";
        type_        =   3;
        mass_        =   ME;
        lifetime_    =   -1;
    }
    else if(StartsWith(name,"mu"))
    {
        name_        =   "mu";
        type_        =   1;
        mass_        =   MMU;
        lifetime_    =   LMU;
    }
    else
    {
        InitByName(name);
    }
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


Particle::Particle(int parent_particle_id,
                   int particle_id,
                   string name,
                   double x,
                   double y,
                   double z,
                   double theta,
                   double phi,
                   double energy,
                   double t,
                   double prop_dist)

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
    ,charge_                ( 1 )
    ,name_                  ( "mu" )
    ,low_                   ( mass_ )
    ,type_                  ( 1 )
    ,parent_particle_id_    ( parent_particle_id )
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
    InitByName(name);
    SetEnergy(energy);
    Location(t, x, y, z, theta, phi);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Particle::Particle(string name,
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
    ,charge_                ( 1 )
    ,name_                  ( name )
    ,low_                   ( mass_ )
    ,type_                  ( 1 )
    ,parent_particle_id_    ( 0 )
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
    if(StartsWith(name,"tau"))
    {
        name_        =   "tau";
        type_        =   2;
        mass_        =   MTAU;
        lifetime_    =   LTAU;
        low_         =   mass_;
    }
    else if(StartsWith(name,"e"))
    {
        name_        =   "e";
        type_        =   3;
        mass_        =   ME;
        lifetime_    =   -1;
        low_         =   mass_;
    }
    else if(StartsWith(name,"mu"))
    {
        name_        =   "mu";
        type_        =   1;
        mass_        =   MMU;
        lifetime_    =   LMU;
        low_         =   mass_;
    }
    else
    {
        InitByName(name);
    }
    SetEnergy(energy);
    Location(t, x, y, z, theta, phi);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Particle::Particle(string name)
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
    ,charge_                ( 1 )
    ,name_                  ( name )
    ,low_                   ( mass_ )
    ,type_                  ( 1 )
    ,parent_particle_id_    ( 0 )
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
    if(StartsWith(name,"tau"))
    {
        name_        =   "tau";
        type_        =   2;
        mass_        =   MTAU;
        lifetime_    =   LTAU;
        low_         =   mass_;
    }
    else if(StartsWith(name,"e"))
    {
        name_        =   "e";
        type_        =   3;
        mass_        =   ME;
        lifetime_    =   -1;
        low_         =   mass_;
    }
    else if(StartsWith(name,"mu"))
    {
        name_        =   "mu";
        type_        =   1;
        mass_        =   MMU;
        lifetime_    =   LMU;
        low_         =   mass_;
    }
    else
    {
        InitByName(name);
    }
    SetEnergy(0);
    Location(0, 0, 0, 0, 0, 0);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Particle& Particle::operator=(const Particle &particle){
    if (this != &particle)
    {
      Particle tmp(particle);
      swap(tmp);
    }
    return *this;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Particle::operator==(const Particle &particle) const
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


bool Particle::operator!=(const Particle &particle) const {
  return !(*this == particle);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


ostream& operator<<(ostream& os, Particle const& particle)
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


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Particle::swap(Particle &particle)
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


void Particle::InitByName(std::string aname){

    string name=aname.length()==0?"?":aname[0]=='a'?aname.substr(1):aname;

    if(name.compare("tau")==0 || name.compare("tau-")==0 || name.compare("tau+")==0)
    {
        if(name.compare("tau+")==0)
        {
            type_    =   -33;
        }
        else
        {
            type_    =   -34;
        }

        mass_       =   MTAU;
        lifetime_   =   LTAU;
    }
    else if(name.compare("mu")==0 || name.compare("mu-")==0 || name.compare("mu+")==0)
    {
        if(name.compare("mu+")==0)
        {
            type_    =   -5;
        }
        else
        {
            type_    =   -6;
        }

        mass_       =   MMU;
        lifetime_   =   LMU;
    }
    else if(StartsWith(name,"stau") )
    {
        if((int)(name.find("stau+"))!=-1)
        {
            type_    =   -9131;
        }
        else
        {
            type_    =   -9132;
        }

        try
        {
            mass_   =   strtod(name.substr(5).c_str(),NULL)*1.e3;
        }
        catch(exception &e)
        {
            mass_   =   0;
        }

        if(mass_<=0)
        {
            mass_   =   MSTAU;
        }

        lifetime_   =   LSTAU;
    }
    else if(name.compare("e")==0 || name.compare("e-")==0 || name.compare("e+")==0)
    {
        if(name.compare("e+")==0)
        {
            type_    =   -2;
        }
        else
        {
            type_    =   -3;
        }

        mass_       =   ME;
        lifetime_   =   -1;
    }
    else if((int)(name.find("nu_"))!=-1)
    {
        if(name.compare("nu_e")==0)
        {
            type_    =   -201;
        }
        else if(name.compare("~nu_e")==0)
        {
            type_    =   -204;
        }
        else if(name.compare("nu_mu")==0)
        {
            type_    =   -202;
        }
        else if(name.compare("~nu_mu")==0)
        {
            type_    =   -205;
        }
        else if(name.compare("nu_tau")==0)
        {
            type_    =   -203;
        }
        else if(name.compare("~nu_tau")==0)
        {
            type_    =   -206;
        }

        mass_       =    0;
        lifetime_   =   -1;
    }
    else
    {
        if(name.compare("delta")==0)
        {
            type_    =   -1002;
        }
        else if(name.compare("brems")==0)
        {
            type_    =   -1001;
        }
        else if(name.compare("munu")==0)
        {
            type_    =   -1004;
        }
        else if(name.compare("epair")==0)
        {
            type_    =   -1003;
        }
        else if(name.compare("hadr")==0)
        {
            type_    =   -1006;
        }
        else if(name.compare("conti")==0)
        {
            type_    =   -1111;
        }
        else
        {
            type_    =   0;
        }

        mass_       =   0;
        lifetime_   =   0;
    }
    name_   =   name;
    low_    =   mass_;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Particle::Location(double time,
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


void Particle::SetEnergy(double e)
{
    energy_             =   e;

    if(energy_ < mass_) energy_ = mass_;

    square_momentum_    =   e*e-mass_*mass_;
    momentum_           =   sqrt(max(square_momentum_,0.0));
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Particle::SetTheta(double theta)
{
    theta_  =   theta;
    theta_ *=   (PI/180);

    costh_  =   cos(theta_);
    sinth_  =   sin(theta_);

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Particle::SetPhi(double phi)
{
    phi_   =    phi;
    phi_  *=    (PI/180);

    cosph_ =    cos(phi_);
    sinph_ =    sin(phi_);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Particle::SetMomentum(double momentum)
{
    momentum_           =   momentum;
    square_momentum_    =   momentum_*momentum_;
    energy_             =   sqrt(square_momentum_ + mass_*mass_);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Particle::SetPropagatedDistance(double prop_dist){
    propagated_distance_ = prop_dist;
}

void Particle::SetX(double x){
    x_ = x;
}

void Particle::SetY(double y){
    y_ = y;
}

void Particle::SetZ(double z){
    z_ = z;
}

void Particle::SetT(double t){
    t_ = t;
}

void Particle::SetMass(double mass){
    mass_ = mass;
}

void Particle::SetLifetime(double lifetime){
    lifetime_ = lifetime;
}

void Particle::SetCharge(double charge){
    charge_ = charge;
}

void Particle::SetName(std::string name){
    name_ = name;
}

void Particle::SetLow(double low){
    low_ = low;
}

void Particle::SetType(int type){
    type_ = type;
}

void Particle::SetParentParticleId(int parent_particle_id){
    parent_particle_id_ = parent_particle_id;
}

void Particle::SetParticleId(int particle_id){
    particle_id_ = particle_id;
}

void Particle::SetXi(double xi){
    xi_ = xi;
}

void Particle::SetYi(double yi){
    yi_ = yi;
}

void Particle::SetZi(double zi){
    zi_ = zi;
}

void Particle::SetTi(double ti){
    ti_ = ti;
}

void Particle::SetEi(double ei){
    ei_ = ei;
}

void Particle::SetXf(double xf){
    xf_ = xf;
}

void Particle::SetYf(double yf){
    yf_ = yf;
}

void Particle::SetZf(double zf){
    zf_ = zf;
}

void Particle::SetTf(double tf){
    tf_ = tf;
}

void Particle::SetEf(double ef){
    ef_ = ef;
}

void Particle::SetXc(double xc){
    xc_ = xc;
}

void Particle::SetYc(double yc){
    yc_ = yc;
}

void Particle::SetZc(double zc){
    zc_ = zc;
}

void Particle::SetTc(double tc){
    tc_ = tc;
}

void Particle::SetEc(double ec){
    ec_ = ec;
}

void Particle::SetElost(double elost){
    elost_ = elost;
}


