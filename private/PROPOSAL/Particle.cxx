/*! \file   Particle.cxx
*   \brief  Source file for the Particle routines.
*
*   For more details see the class documentation.
*
*   \date   21.06.2010
*   \author Jan-Hendrik Köhne
*/


#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include "PROPOSAL/Particle.h"
#include "PROPOSAL/methods.h"



using namespace std;

Particle::Particle()
:r  (0)
,x  (0)
,y  (0)
,z  (0)
,t  (0)
,theta  (0)
,phi    (0)
,costh  (1.)
,sinth  (0)
,cosph  (1.)
,sinph  (0.)
,p  (0)
,p2 (0)
,e  (0)
,m  (MMU)
,l  (LMU)
,c  (1)
,name   ("mu")
,low    (m)
,type   (1)
,igen   (0)
,gens    (1)
,xi     (0)
,yi     (0)
,zi     (0)
,ti     (0)
,Ei     (0)
,xf     (0)
,yf     (0)
,zf     (0)
,tf     (0)
,Ef     (0)
,xc     (0)
,yc     (0)
,zc     (0)
,tc     (0)
,Ec     (0)
,Elost  (0)
,df (false)
,jt (false)
{
}

//----------------------------------------------------------------------------//



Particle::Particle(int Igen,
                   int Gens,
                   string Name,
                   double X,
                   double Y,
                   double Z,
                   double Theta,
                   double Phi,
                   double E,
                   double T,
                   double R,
                   Particle *P)
:r  (R)
,x  (0)
,y  (0)
,z  (0)
,t  (0)
,theta  (0)
,phi    (0)
,costh  (1.)
,sinth  (0)
,cosph  (1.)
,sinph  (0.)
,p  (0)
,p2 (0)
,e  (0)
,m  (MMU)
,l  (LMU)
,c  (1)
,name   ("mu")
,low    (m)
,type   (1)
,igen   (Igen)
,gens   (Gens)
,xi     (0)
,yi     (0)
,zi     (0)
,ti     (0)
,Ei     (0)
,xf     (0)
,yf     (0)
,zf     (0)
,tf     (0)
,Ef     (0)
,xc     (0)
,yc     (0)
,zc     (0)
,tc     (0)
,Ec     (0)
,Elost  (0)
,df (false)
,jt (false)
{
    initByName(Name);
    setEnergy(E);
    location(Name, T, X, Y, Z, Theta, Phi);


    if(P!=NULL)
    {
        xi      =   P->xi;
        yi      =   P->yi;
        zi      =   P->zi;
        ti      =   P->ti;
        Ei      =   P->Ei;
        xf      =   P->xf;
        yf      =   P->yf;
        zf      =   P->zf;
        tf      =   P->tf;
        Ef      =   P->Ef;
        xc      =   P->xc;
        yc      =   P->yc;
        zc      =   P->zc;
        tc      =   P->tc;
        Ec      =   P->Ec;
        Elost   =   P->Elost;
    }
}


//----------------------------------------------------------------------------//


Particle::Particle(int Igen,
                   int Gens,
                   string Name,
                   double X,
                   double Y,
                   double Z,
                   double Theta,
                   double Phi,
                   double E,
                   double T,
                   double R)
    :r  (R)
    ,x  (0)
    ,y  (0)
    ,z  (0)
    ,t  (0)
    ,theta  (0)
    ,phi    (0)
    ,costh  (1.)
    ,sinth  (0)
    ,cosph  (1.)
    ,sinph  (0.)
    ,p  (0)
    ,p2 (0)
    ,e  (0)
    ,m  (MMU)
    ,l  (LMU)
    ,c  (1)
    ,name   ("mu")
    ,low    (m)
    ,type   (1)
    ,igen   (Igen)
    ,gens   (Gens)
    ,xi     (0)
    ,yi     (0)
    ,zi     (0)
    ,ti     (0)
    ,Ei     (0)
    ,xf     (0)
    ,yf     (0)
    ,zf     (0)
    ,tf     (0)
    ,Ef     (0)
    ,xc     (0)
    ,yc     (0)
    ,zc     (0)
    ,tc     (0)
    ,Ec     (0)
    ,Elost  (0)
    ,df (false)
    ,jt (false)
{
    initByName(Name);
    setEnergy(E);
    location(Name, T, X, Y, Z, Theta, Phi);
}


//----------------------------------------------------------------------------------------------------//


Particle::Particle(string aname,
                   double X,
                   double Y,
                   double Z,
                   double Theta,
                   double Phi,
                   double E,
                   double T)
:r  (0)
,x  (0)
,y  (0)
,z  (0)
,t  (0)
,theta  (0)
,phi    (0)
,costh  (1.)
,sinth  (0)
,cosph  (1.)
,sinph  (0.)
,p  (0)
,p2 (0)
,e  (0)
,m  (MMU)
,l  (LMU)
,c  (1)
,name   ("mu")
,low    (m)
,type   (1)
,igen   (0)
,gens   (1)
,xi     (0)
,yi     (0)
,zi     (0)
,ti     (0)
,Ei     (0)
,xf     (0)
,yf     (0)
,zf     (0)
,tf     (0)
,Ef     (0)
,xc     (0)
,yc     (0)
,zc     (0)
,tc     (0)
,Ec     (0)
,Elost  (0)
,df (false)
,jt (false)
{
    initByName(aname);
    setEnergy(E);
    location(aname, T, X, Y, Z, Theta, Phi);
}

//----------------------------------------------------------------------------//

void Particle::initByName(std::string aname){
    string name=aname.length()==0?"?":aname[0]=='a'?aname.substr(1):aname;

    if(name.compare("tau")==0 || name.compare("tau-")==0 || name.compare("tau+")==0)
    {
        if(name.compare("tau+")==0)
        {
            type    =   -33;
        }
        else
        {
            type    =   -34;
        }

        m   =   MTAU;
        l   =   LTAU;
    }
    else if(name.compare("mu")==0 || name.compare("mu-")==0 || name.compare("mu+")==0)
    {
        if(name.compare("mu+")==0)
        {
            type    =   -5;
        }
        else
        {
            type    =   -6;
        }

        m   =   MMU;
        l   =   LMU;
    }
    else if(StartsWith(name,"stau") )
    {
        if((int)(name.find("stau+"))!=-1)
        {
            type    =   -9131;
        }
        else
        {
            type    =   -9132;
        }

        try
        {
            m   =   strtod(name.substr(5).c_str(),NULL)*1.e3;
        }
        catch(exception &e)
        {
            m   =   0;
        }

        if(m<=0)
        {
            m   =   MSTAU;
        }

        l   =   LSTAU;
    }
    else if(name.compare("e")==0 || name.compare("e-")==0 || name.compare("e+")==0)
    {
        if(name.compare("e+")==0)
        {
            type    =   -2;
        }
        else
        {
            type    =   -3;
        }

        m   =   ME;
        l   =   -1;
    }
    else if((int)(name.find("nu_"))!=-1)
    {
        if(name.compare("nu_e")==0)
        {
            type    =   -201;
        }
        else if(name.compare("~nu_e")==0)
        {
            type    =   -204;
        }
        else if(name.compare("nu_mu")==0)
        {
            type    =   -202;
        }
        else if(name.compare("~nu_mu")==0)
        {
            type    =   -205;
        }
        else if(name.compare("nu_tau")==0)
        {
            type    =   -203;
        }
        else if(name.compare("~nu_tau")==0)
        {
            type    =   -206;
        }

        m   =   0;
        l   =   -1;
    }
    else
    {
        if(name.compare("delta")==0)
        {
            type    =   -1002;
        }
        else if(name.compare("brems")==0)
        {
            type    =   -1001;
        }
        else if(name.compare("munu")==0)
        {
            type    =   -1004;
        }
        else if(name.compare("epair")==0)
        {
            type    =   -1003;
        }
        else if(name.compare("hadr")==0)
        {
            type    =   -1006;
        }
        else if(name.compare("conti")==0)
        {
            type    =   -1111;
        }
        else
        {
            type    =   0;
        }

        m   =   0;
        l   =   0;
    }

    low =   m;
}

//----------------------------------------------------------------------------//

void Particle::location(string name,
                        double time,
                        double x,
                        double y,
                        double z,
                        double theta,
                        double phi)
{
    this->name  =   name;
    r           =   0;
    t           =   time;
    this->x     =   x;
    this->y     =   y;
    this->z     =   z;
    this->theta =   theta;
    this->phi   =   phi;
    theta       *=  (PI/180);
    phi         *=  (PI/180);
    costh       =   cos(theta);
    sinth       =   sin(theta);
    cosph       =   cos(phi);
    sinph       =   sin(phi);
}

//----------------------------------------------------------------------------//

//Setter

void Particle::SetEnergy(double e)
{
    energy_             =   e;
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
    theta_ = theta;
    theta_ *=  (PI/180);

    costh_  =   cos(theta);
    sinth_  =   sin(theta);

}
//----------------------------------------------------------------------------//
void Particle::SetPhi(double phi)
{
    phi_ = phi;
    phi  *=  (PI/180);

    cosph =   cos(phi);
    sinph =   sin(phi);
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


