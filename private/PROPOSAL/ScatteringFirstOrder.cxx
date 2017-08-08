/*! \file   ScatteringFirstOrder.cxx
*   \brief  Source file for the ScatteringFirstOrder routines.
*
*   For more details see the class documentation.
*
*   \date   2013.06.13
*   \author Tomasz Fuchs
*/


// #include <cmath>
// #include <algorithm>
// #include <stdlib.h>

#include "PROPOSAL/methods.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/ScatteringFirstOrder.h"
// #include "PROPOSAL/Output.h"

using namespace std;
using namespace PROPOSAL;


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------public member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void ScatteringFirstOrder::Scatter(double dr, PROPOSALParticle* part, Medium* med)
{
    double theta0, rnd1, rnd2, sx, tx, sy, ty, sz, tz;

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


    Vector3D position;
    Vector3D direction;

    long double sinth, costh,sinph,cosph;
    sinth = (long double) sin(part->GetDirection().GetTheta());
    costh = (long double) cos(part->GetDirection().GetTheta());
    sinph = (long double) sin(part->GetDirection().GetPhi());
    cosph = (long double) cos(part->GetDirection().GetPhi());

    position = part->GetPosition();

    // Rotation towards all tree axes
    direction = sz*part->GetDirection();
    direction = direction + sx*Vector3D(costh*cosph, costh*sinph, -sinth);
    direction = direction + sy*Vector3D(-sinph, cosph, 0.);

    position = position + dr*direction;

    // Rotation towards all tree axes
    direction = tz*part->GetDirection();
    direction = direction + tx*Vector3D(costh*cosph, costh*sinph, -sinth);
    direction = direction + ty*Vector3D(-sinph, cosph, 0.);

    direction.CalculateSphericalCoordinates();

    part->SetPosition(position);
    part->SetDirection(direction);

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


double ScatteringFirstOrder::CalculateTheta0(double dr, PROPOSALParticle* part, Medium* med)
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





