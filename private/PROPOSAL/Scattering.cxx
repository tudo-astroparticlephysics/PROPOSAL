/*! \file   Scattering.cxx
*   \brief  Source file for the scattering routines.
*
*   For more details see the class documentation.
*
*   \date   2013.06.13
*   \author Tomasz Fuchs
*/


#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Scattering.h"
#include "PROPOSAL/Output.h"

using namespace std;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Scattering::Scattering( )
{
    standard_normal_    =   new StandardNormal(IROMB, IMAXS, IPREC);
}


Scattering::Scattering(StandardNormal* standard_normal)
{
    if(standard_normal != NULL)
    {
        standard_normal_    =   standard_normal;
    }
    else
    {
        standard_normal_    =   new StandardNormal(IROMB, IMAXS, IPREC);
    }
}

Scattering::Scattering(const Scattering &scattering)
{
    if(scattering.standard_normal_ != NULL)
    {
        standard_normal_ = new StandardNormal( *scattering.standard_normal_);
    }
    else
    {
        standard_normal_    =   new StandardNormal(IROMB, IMAXS, IPREC);
    }
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Scattering& Scattering::operator=(const Scattering &scattering){
    if (this != &scattering)
    {
      Scattering tmp(scattering);
      swap(tmp);
    }
    return *this;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Scattering::operator==(const Scattering &scattering) const
{
    if( standard_normal_ != NULL && scattering.standard_normal_ != NULL)
    {
        if( *standard_normal_   != *scattering.standard_normal_)                                        return false;
    }
    return true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Scattering::operator!=(const Scattering &scattering) const {
  return !(*this == scattering);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Scattering::swap(Scattering &scattering)
{
    using std::swap;

    if(scattering.standard_normal_ != NULL)
    {
        standard_normal_->swap(*scattering.standard_normal_) ;
    }
    else
    {
        standard_normal_ = NULL;
    }
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//------------------------private member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Scattering::CalculateTheta0(double dr, Particle* part, Medium* med)
{
    double y = dr/med->GetRadiationLength();
    double beta = 1./sqrt(1 +  part->GetMass() * part->GetMass()/ (part->GetMomentum()*part->GetMomentum() ));
    y = 13.6/(part->GetMomentum()* beta ) *sqrt(y)*( 1.+0.038*log(y) );
    return y;
}


void Scattering::Scatter(double dr, Particle* part, Medium* med)
{
        double theta0, theta_max, rnd1, rnd2, sx, tx, sy, ty, sz, tz, ax, ay, az;
        double x,y,z;

        theta0     =   CalculateTheta0(dr, part,med);

        theta_max    =   1./SQRT2;


        rnd1    =   standard_normal_-> StandardNormalRandomNumber(RandomDouble(), 0, theta0, -theta_max, theta_max, false);
        rnd2    =   standard_normal_-> StandardNormalRandomNumber(RandomDouble(), 0, theta0, -theta_max, theta_max, false);
        sx      =   (rnd1/SQRT3+rnd2)/2;
        tx      =   rnd2;

        rnd1    =   standard_normal_-> StandardNormalRandomNumber(RandomDouble(), 0, theta0, -theta_max, theta_max, false);
        rnd2    =   standard_normal_-> StandardNormalRandomNumber(RandomDouble(), 0, theta0, -theta_max, theta_max, false);

        sy      =   (rnd1/SQRT3+rnd2)/2;
        ty      =   rnd2;

        sz      =   sqrt(max(1.-(sx*sx+sy*sy), 0.));
        tz      =   sqrt(max(1.-(tx*tx+ty*ty), 0.));

        double sinth, costh,sinph,cosph;
        double theta, phi;

        sinth = part->GetSinTheta();
        costh = part->GetCosTheta();
        sinph = part->GetSinPhi();
        cosph = part->GetCosPhi();

        x   = part->GetX();
        y   = part->GetY();
        z   = part->GetZ();


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
        sinth   =   sqrt(max(1-costh*costh, 0.));

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

        part->SetX(x);
        part->SetY(y);
        part->SetZ(z);
        part->SetPhi(phi);
        part->SetTheta(theta);
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//





