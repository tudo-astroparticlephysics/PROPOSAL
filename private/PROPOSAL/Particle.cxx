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


using namespace std;

Particle::Particle()
    :propagation_distance_  (0)
    ,x_                     (0)
    ,y_                     (0)
    ,z_                     (0)
    ,t_                     (0)
    ,theta_                 (0)
    ,phi_                   (0)
    ,costh_                 (1.)
    ,sinth_                 (0)
    ,cosph_                 (1.)
    ,sinph_                 (0.)
    ,momentum_              (0)
    ,square_momentum_       (0)
    ,energy_                (0)
    ,mass_                  (MMU)
    ,lifetime_              (LMU)
    ,charge_                (1)
    ,name_                  ("mu")
    ,low_                   (mass_)
    ,type_                  (1)
    ,parent_particle_id_    (0)
    ,particle_id_           (1)
    ,xi_                    (0)
    ,yi_                    (0)
    ,zi_                    (0)
    ,ti_                    (0)
    ,ei_                    (0)
    ,xf_                    (0)
    ,yf_                    (0)
    ,zf_                    (0)
    ,tf_                    (0)
    ,ef_                    (0)
    ,xc_                    (0)
    ,yc_                    (0)
    ,zc_                    (0)
    ,tc_                    (0)
    ,ec_                    (0)
    ,elost_                 (0)
{
}

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

    :propagation_distance_  (prop_dist)
    ,x_                     (x)
    ,y_                     (y)
    ,z_                     (z)
    ,t_                     (t)
    ,theta_                 (theta)
    ,phi_                   (phi)
    ,momentum_              (0)
    ,square_momentum_       (0)
    ,energy_                (0)
    ,mass_                  (MMU)
    ,lifetime_              (LMU)
    ,charge_                (1)
    ,name_                  ("mu")
    ,low_                   (mass_)
    ,type_                  (1)
    ,parent_particle_id_    (parent_particle_id)
    ,particle_id_           (particle_id)
    ,xi_                    (0)
    ,yi_                    (0)
    ,zi_                    (0)
    ,ti_                    (0)
    ,ei_                    (0)
    ,xf_                    (0)
    ,yf_                    (0)
    ,zf_                    (0)
    ,tf_                    (0)
    ,ef_                    (0)
    ,xc_                    (0)
    ,yc_                    (0)
    ,zc_                    (0)
    ,tc_                    (0)
    ,ec_                    (0)
    ,elost_                 (0)
{
    InitByName(name);
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

    :propagation_distance_  (prop_dist)
    ,x_                     (x)
    ,y_                     (y)
    ,z_                     (z)
    ,t_                     (t)
    ,theta_                 (theta)
    ,phi_                   (phi)
    ,momentum_              (0)
    ,square_momentum_       (0)
    ,energy_                (energy)
    ,mass_                  (MMU)
    ,lifetime_              (LMU)
    ,charge_                (1)
    ,name_                  ("mu")
    ,low_                   (mass_)
    ,type_                  (1)
    ,parent_particle_id_    (parent_particle_id)
    ,particle_id_           (particle_id)
    ,xi_                    (0)
    ,yi_                    (0)
    ,zi_                    (0)
    ,ti_                    (0)
    ,ei_                    (0)
    ,xf_                    (0)
    ,yf_                    (0)
    ,zf_                    (0)
    ,tf_                    (0)
    ,ef_                    (0)
    ,xc_                    (0)
    ,yc_                    (0)
    ,zc_                    (0)
    ,tc_                    (0)
    ,ec_                    (0)
    ,elost_                 (0)

{
    InitByName(name);
    SetEnergy(energy);
    Location(t, x, y, z, theta, phi);
}


//----------------------------------------------------------------------------------------------------//


Particle::Particle(string name,
                   double x,
                   double y,
                   double z,
                   double theta,
                   double phi,
                   double energy,
                   double t)

    :propagation_distance_  (0)
    ,x_                     (x)
    ,y_                     (y)
    ,z_                     (z)
    ,t_                     (t)
    ,theta_                 (theta)
    ,phi_                   (phi)
    ,momentum_              (0)
    ,square_momentum_       (0)
    ,energy_                (energy)
    ,mass_                  (MMU)
    ,lifetime_              (LMU)
    ,charge_                (1)
    ,name_                  (name)
    ,low_                   (mass_)
    ,type_                  (1)
    ,parent_particle_id_    (0)
    ,particle_id_           (1)
    ,xi_                    (0)
    ,yi_                    (0)
    ,zi_                    (0)
    ,ti_                    (0)
    ,ei_                    (0)
    ,xf_                    (0)
    ,yf_                    (0)
    ,zf_                    (0)
    ,tf_                    (0)
    ,ef_                    (0)
    ,xc_                    (0)
    ,yc_                    (0)
    ,zc_                    (0)
    ,tc_                    (0)
    ,ec_                    (0)
    ,elost_                 (0)
{
    InitByName(name);
    SetEnergy(energy);
    Location(t, x, y, z, theta, phi);
}

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

void Particle::Location(double time,
                        double x,
                        double y,
                        double z,
                        double theta,
                        double phi)
{
    //propagation_distance_ = 0; <--- seems to be wrong
    t_           =   time;
    x_           =   x;
    y_           =   y;
    z_           =   z;
    theta_       =   theta;
    phi_         =   phi;
    theta_       *=  (PI/180);   //<--deg or rad??
    phi_         *=  (PI/180);   //<--deg or rad??
    costh_       =   cos(theta);
    sinth_       =   sin(theta);
    cosph_       =   cos(phi);
    sinph_       =   sin(phi);
}

//----------------------------------------------------------------------------//

//Setter

void Particle::SetEnergy(double e)
{
    energy_             =   e;

    if(energy_ < mass_) energy_ = mass_;

    square_momentum_    =   e*e-mass_*mass_;
    momentum_           =   sqrt(max(square_momentum_,0.0));
}

//----------------------------------------------------------------------------//
void Particle::SetPropagationDistance(double prop_dist)
{
    propagation_distance_ = prop_dist;
}
//----------------------------------------------------------------------------//
void Particle::SetX(double x)
{
    x_ = x;
}
//----------------------------------------------------------------------------//
void Particle::SetY(double y)
{
    y_ = y;
}
//----------------------------------------------------------------------------//
void Particle::SetZ(double z)
{
    z_ = z;
}
//----------------------------------------------------------------------------//
void Particle::SetT(double t)
{
    t_ = t;
}
//----------------------------------------------------------------------------//
void Particle::SetTheta(double theta)
{
    theta_  =   theta;
    theta_ *=   (PI/180);

    costh_  =   cos(theta);
    sinth_  =   sin(theta);

}
//----------------------------------------------------------------------------//
void Particle::SetPhi(double phi)
{
    phi_   =    phi;
    phi_  *=    (PI/180);

    cosph_ =    cos(phi);
    sinph_ =    sin(phi);
}
//----------------------------------------------------------------------------//
void Particle::SetMomentum(double momentum)
{
    momentum_           =   momentum;
    square_momentum_    =   momentum_*momentum_;
    energy_             =   sqrt(square_momentum_ + mass_*mass_);
}
//----------------------------------------------------------------------------//
void Particle::SetMass(double mass)
{
    mass_ = mass;
}
//----------------------------------------------------------------------------//
void Particle::SetLifetime(double lifetime)
{
    lifetime_ = lifetime;
}
//----------------------------------------------------------------------------//
void Particle::SetCharge(double charge)
{
    charge_ = charge;
}
//----------------------------------------------------------------------------//
void Particle::SetName(std::string name)
{
    name_ = name;
}
//----------------------------------------------------------------------------//
void Particle::SetLow(double low)
{
    low_ = low;
}
//----------------------------------------------------------------------------//
void Particle::SetType(int type)
{
    type_ = type;
}
//----------------------------------------------------------------------------//
void Particle::SetParentParticleId(int parent_particle_id)
{
    parent_particle_id_ = parent_particle_id;
}
//----------------------------------------------------------------------------//
void Particle::SetParticleId(int particle_id)
{
    particle_id_ = particle_id;
}
//----------------------------------------------------------------------------//
void Particle::SetXi(double xi)
{
    xi_ = xi;
}
//----------------------------------------------------------------------------//
void Particle::SetYi(double yi)
{
    yi_ = yi;
}
//----------------------------------------------------------------------------//
void Particle::SetZi(double zi)
{
    zi_ = zi;
}
//----------------------------------------------------------------------------//
void Particle::SetTi(double ti)
{
    ti_ = ti;
}
//----------------------------------------------------------------------------//
void Particle::SetEi(double ei)
{
    ei_ = ei;
}
//----------------------------------------------------------------------------//
void Particle::SetXf(double xf)
{
    xf_ = xf;
}
//----------------------------------------------------------------------------//
void Particle::SetYf(double yf)
{
    yf_ = yf;
}
//----------------------------------------------------------------------------//
void Particle::SetZf(double zf)
{
    zf_ = zf;
}
//----------------------------------------------------------------------------//
void Particle::SetTf(double tf)
{
    tf_ = tf;
}
//----------------------------------------------------------------------------//
void Particle::SetEf(double ef)
{
    ef_ = ef;
}
//----------------------------------------------------------------------------//
void Particle::SetXc(double xc)
{
    xc_ = xc;
}
//----------------------------------------------------------------------------//
void Particle::SetYc(double yc)
{
    yc_ = yc;
}
//----------------------------------------------------------------------------//
void Particle::SetZc(double zc)
{
    zc_ = zc;
}
//----------------------------------------------------------------------------//
void Particle::SetTc(double tc)
{
    tc_ = tc;
}
//----------------------------------------------------------------------------//
void Particle::SetEc(double ec)
{
    ec_ = ec;
}
//----------------------------------------------------------------------------//
void Particle::SetElost(double elost)
{
    elost_ = elost;
}


