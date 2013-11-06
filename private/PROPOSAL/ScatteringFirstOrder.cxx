/*! \file   ScatteringFirstOrder.cxx
*   \brief  Source file for the ScatteringFirstOrder routines.
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
#include "PROPOSAL/ScatteringFirstOrder.h"
#include "PROPOSAL/Output.h"

#include <boost/math/special_functions/erf.hpp>
#define erfInv(x)   boost::math::erf_inv(x)

using namespace std;


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------public member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void ScatteringFirstOrder::Scatter(double dr, Particle* part, Medium* med)
{
        double theta0, rnd1, rnd2, sx, tx, sy, ty, sz, tz, ax, ay, az;
        double x,y,z;

        theta0     =   CalculateTheta0(dr, part,med);

        rnd1 = SQRT2*theta0*erfInv( 2.*(RandomDouble()-0.5) );
        rnd2 = SQRT2*theta0*erfInv( 2.*(RandomDouble()-0.5) );

        sx      =   (rnd1/SQRT3+rnd2)/2;
        tx      =   rnd2;

        rnd1 = SQRT2*theta0*erfInv(2*(RandomDouble()-0.5));
        rnd2 = SQRT2*theta0*erfInv(2*(RandomDouble()-0.5));

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
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

ScatteringFirstOrder::ScatteringFirstOrder( )
{

}



//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//------------------------private member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double ScatteringFirstOrder::CalculateTheta0(double dr, Particle* part, Medium* med)
{
    double y = dr/med->GetRadiationLength();
    double beta = 1./sqrt(1 +  part->GetMass() * part->GetMass()/ (part->GetMomentum()*part->GetMomentum() ));
    y = 13.6/(part->GetMomentum()* beta ) *sqrt(y)*( 1.+0.088*log10(y) );
    return y;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//





