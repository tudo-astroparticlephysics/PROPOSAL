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
#include "PROPOSAL/math/MathModel.h"
#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/scattering/ScatteringFirstOrder.h"
#include "PROPOSAL/Output.h"

using namespace std;
using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //


ScatteringFirstOrder::ScatteringFirstOrder(PROPOSALParticle& particle, const Medium& medium)
    : Scattering(particle)
    , medium_(medium.clone())
{
}

ScatteringFirstOrder::ScatteringFirstOrder(const ScatteringFirstOrder& scattering)
    : Scattering(scattering)
    , medium_(scattering.medium_->clone())
{
}

ScatteringFirstOrder::ScatteringFirstOrder(PROPOSALParticle& particle, const ScatteringFirstOrder& scattering)
    : Scattering(particle)
    , medium_(scattering.medium_->clone())
{
}

ScatteringFirstOrder::~ScatteringFirstOrder()
{
    delete medium_;
}

// ------------------------------------------------------------------------- //
// Private methods
// ------------------------------------------------------------------------- //

double ScatteringFirstOrder::CalculateTheta0(double dr)
{
    double y = dr/medium_->GetRadiationLength();
    double beta = 1./sqrt(1 +  particle_.GetMass() * particle_.GetMass()/ (particle_.GetMomentum()*particle_.GetMomentum() ));
    y = 13.6/(particle_.GetMomentum()* beta ) *sqrt(y)*( 1.+0.088*log10(y) );
    return y;
}

//----------------------------------------------------------------------------//

Scattering::RandomAngles ScatteringFirstOrder::CalculateRandomAngle(double dr, double ei, double ef)
{
    (void)ei;
    (void)ef;

    double Theta0,rnd1,rnd2;
    Scattering::RandomAngles random_angles;

    Theta0 = CalculateTheta0(dr);

    rnd1 = SQRT2*Theta0*erfInv( 2.*(RandomGenerator::Get().RandomDouble()-0.5) );
    rnd2 = SQRT2*Theta0*erfInv( 2.*(RandomGenerator::Get().RandomDouble()-0.5) );

    random_angles.sx = (rnd1/SQRT3+rnd2)/2;
    random_angles.tx = rnd2;

    rnd1 = SQRT2*Theta0*erfInv(2*(RandomGenerator::Get().RandomDouble()-0.5));
    rnd2 = SQRT2*Theta0*erfInv(2*(RandomGenerator::Get().RandomDouble()-0.5));

    random_angles.sy = (rnd1/SQRT3+rnd2)/2;
    random_angles.ty = rnd2;

    return random_angles;
}
