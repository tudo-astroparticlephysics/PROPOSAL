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


ScatteringFirstOrder::ScatteringFirstOrder()
{
}

ScatteringFirstOrder::ScatteringFirstOrder(const ScatteringFirstOrder& scattering)
    : Scattering(scattering)
{
}

ScatteringFirstOrder::~ScatteringFirstOrder()
{
}

// ------------------------------------------------------------------------- //
// Public methods
// ------------------------------------------------------------------------- //


void ScatteringFirstOrder::EnableInterpolation(const PROPOSALParticle& particle,
                                            const std::vector<CrossSection*>& cross_sections,
                                            std::string filepath)
{
    (void)particle;
    (void)cross_sections;
    (void)filepath;

    log_warn("No interpolation implemented for ScatteringFirstOrder");
}

void ScatteringFirstOrder::DisableInterpolation()
{
    log_warn("No interpolation implemented for ScatteringFirstOrder");
}

// ------------------------------------------------------------------------- //
// Private methods
// ------------------------------------------------------------------------- //

double ScatteringFirstOrder::CalculateTheta0(const PROPOSALParticle& particle, const Medium& med, double dr)
{
    double y = dr/med.GetRadiationLength();
    double beta = 1./sqrt(1 +  particle.GetMass() * particle.GetMass()/ (particle.GetMomentum()*particle.GetMomentum() ));
    y = 13.6/(particle.GetMomentum()* beta ) *sqrt(y)*( 1.+0.088*log10(y) );
    return y;
}

//----------------------------------------------------------------------------//

Scattering::RandomAngles ScatteringFirstOrder::CalculateRandomAngle(const PROPOSALParticle& particle, const std::vector<CrossSection*>& cross_sections, double dr, double ei, double ef)
{
    (void)ei;
    (void)ef;

    double Theta0,rnd1,rnd2;
    Scattering::RandomAngles random_angles;

    Theta0 = CalculateTheta0(particle, cross_sections.at(0)->GetParametrization().GetMedium(), dr);

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
