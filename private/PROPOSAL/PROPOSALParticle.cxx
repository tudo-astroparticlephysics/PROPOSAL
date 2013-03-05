/*! \file   PROPOSALParticle.cxx
*   \brief  Source file for the PROPOSALParticle routines.
*
*   For more details see the class documentation.
*
*   \date   21.06.2010
*   \author Jan-Hendrik Köhne
*/


#include "PROPOSAL/methods.h"
#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include <exception>

#include "PROPOSAL/PROPOSALParticle.h"

using namespace std;


PROPOSALParticle::PROPOSALParticle()
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


PROPOSALParticle::PROPOSALParticle(Propagate *pr, string name)
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
    if(StartsWith(name,"tau"))
    {
        this->name  =   "tau";
        type        =   2;
        m           =   MTAU;
        l           =   LTAU;
    }
    else if(StartsWith(name,"e"))
    {
        this->name  =   "e";
        type        =   3;
        m           =   ME;
        l           =   -1;
    }
    else if(StartsWith(name,"monopole"))
    {
        type        =   4;

        try
        {
            m=strtod(name.substr(9).c_str(),NULL)*1.e3;
        }
        catch(exception &e)
        {
            m       =   0;
        }

        if(m<=0)
        {
            m       =   MMON;
        }

        this->name  =   "monopole-"+output->f(m*1.e-3);
        c           =   CMON;
        l           =   -1;
    }
    else if(StartsWith(name,"stau"))
    {
        type        =   5;

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

        this->name  =   "stau-"+output->f(m*1.e-3);
        l           =   LSTAU;
    }
    else
    {
        this->name  =   "mu";
        type        =   1;
        m           =   MMU;
        l           =   LMU;
    }

    low                 =   m;
    e                   =   m;
    this->propagate_    =   pr;

    if(PhysicsModel::get_elow()>low)
    {
        low=PhysicsModel::get_elow();
    }

    if(PhysicsModel::get_ebig()<100*low)
    {
        PhysicsModel::set_ebig(max(100*low, BIGENERGY));
    }

    scattering_ =   new Scattering(this);
    integral_   =   new Integral(IROMB, IMAXS, IPREC2);
}

//----------------------------------------------------------------------------//



PROPOSALParticle::PROPOSALParticle(int Igen, int Gens, string Name, double X, double Y, double Z, double Theta, double Phi, double E, double T, double R, PROPOSALParticle *P)
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


PROPOSALParticle::PROPOSALParticle(int Igen, int Gens, string Name, double X, double Y, double Z, double Theta, double Phi, double E, double T, double R)
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


PROPOSALParticle::PROPOSALParticle(string aname, double X, double Y, double Z, double Theta, double Phi, double E, double T)
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

void PROPOSALParticle::initByName(std::string aname){
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

void PROPOSALParticle::location(string name, double time, double x, double y, double z, double theta, double phi)
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


void PROPOSALParticle::advance(double dr, double ei, double ef)
{
    long double ax=0, ay=0, az=0;
    long double tho=0, max_=0, rnd1=0, rnd2=0, sx=0, sy=0, sz=0, tx=0, ty=0, tz=0;

    r   +=  dr;

    if(propagate_->get_exactTime())
    {
        t   +=  getdt(ei, ef)/propagate_->get_rho();
    }
    else
    {
        t   +=  dr/SPEED;
    }


    if(propagate_->get_molieScat())
    {

        tho     =   (long double)scattering_->gettho(dr, ei, ef);

        max_    =   1/SQRT2;
        rnd1    =   (long double)propagate_->get_Standard()->sndrn(RandomDouble(), 0, tho, -max_, max_, false);
        rnd2    =   (long double)propagate_->get_Standard()->sndrn(RandomDouble(), 0, tho, -max_, max_, false);

        sx      =   (rnd1/SQRT3+rnd2)/2;
        tx      =   rnd2;
        rnd1    =   (long double)propagate_->get_Standard()->sndrn(RandomDouble(), 0, tho, -max_, max_, false);
        rnd2    =   (long double)propagate_->get_Standard()->sndrn(RandomDouble(), 0, tho, -max_, max_, false);


        sy      =   (rnd1/SQRT3+rnd2)/2;
        ty      =   rnd2;


        sz      =   sqrt(max(1.-(sx*sx+sy*sy), (long double)0.));

        tz      =   sqrt(max(1.-(tx*tx+ty*ty), (long double)0.));


        ax      =   sinth*cosph*sz+costh*cosph*sx-sinph*sy;
        ay      =   sinth*sinph*sz+costh*sinph*sx+cosph*sy;
        az      =   costh*sz-sinth*sx;

        x       +=  ax*dr;
        y       +=  ay*dr;
        z       +=  az*dr;


        ax      =   sinth*cosph*tz+costh*cosph*tx-sinph*ty;
        ay      =   sinth*sinph*tz+costh*sinph*tx+cosph*ty;
        az      =   costh*tz-sinth*tx;



        costh   =   az;
        sinth   =   sqrt(max(1-costh*costh, (long double)0));

        if(sinth!=0)
        {
            sinph   =   ay/sinth;
            cosph   =   ax/sinth;
        }

        if(costh>1)
        {
            theta   =   acos(1)*180./PI;
        }
        else if(costh<-1)
        {
            theta   =   acos(-1)*180./PI;
        }
        else
        {
            theta   =   acos(costh)*180./PI;
        }

        if(cosph>1)
        {
            phi =   acos(1)*180./PI;
        }
        else if(cosph<-1)
        {
            phi =   acos(-1)*180./PI;
        }
        else
        {
            phi =   acos(cosph)*180./PI;
        }

        if(sinph<0)
        {
            phi =   360.-phi;
        }

        if(phi>=360)
        {
            phi -=  360.;
        }

    }
    else
    {
        x   +=  sinth*cosph*dr;
        y   +=  sinth*sinph*dr;
        z   +=  costh*dr;
    }

}

//----------------------------------------------------------------------------//



void PROPOSALParticle::setEnergy(double e)
{
    this->e =   e;
    p2      =   e*e-m*m;
    p       =   sqrt(max(p2,0.0));
}

//----------------------------------------------------------------------------//


double PROPOSALParticle::function(double E)
{
    double aux;

    aux     =   propagate_->get_cros()->function(E);
    aux     *=  e/(p*SPEED);

    return aux;
}

//----------------------------------------------------------------------------//


double PROPOSALParticle::getdt(double ei, double ef)
{
    if(jt)
    {
        if(abs(ei-ef)>abs(ei)*HALF_PRECISION)
        {
            double aux  =   interpolateJ_->interpolate(ei);
            double aux2 =   aux - interpolateJ_->interpolate(ef);

            if(abs(aux2)>abs(aux)*HALF_PRECISION)
            {
                return aux2;
            }
        }

        return interpolateJdf_->interpolate((ei+ef)/2)*(ef-ei);
    }
    else
    {
        return integral_->integrateWithLog(ei, ef, this);
    }
}

//----------------------------------------------------------------------------//


double PROPOSALParticle::functionInt(double e)
{
    if(df)
    {
        return function(e);
    }
    else
    {
        return integral_->integrateWithLog(e, low, this);
    }
}

//----------------------------------------------------------------------------//

double PROPOSALParticle::alpha(double v)
{
    double alpha0   =   0.007297352533285885;
    double Q        =   v*e;

    return alpha0 /( 1- 2*alpha0/(3*PI)*log(Q*Q/ME));

}

