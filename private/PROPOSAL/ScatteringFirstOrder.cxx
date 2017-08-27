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
#include "PROPOSAL/MathModel.h"
#include "PROPOSAL/CrossSections.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/ScatteringFirstOrder.h"
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

void ScatteringFirstOrder::Scatter(PROPOSALParticle& particle,
                                const std::vector<CrossSections*>& cross_sections,
                                double dr,
                                double ei,
                                double ef)
{
    (void) ei;
    (void) ef;

    double theta0, rnd1, rnd2, sx, tx, sy, ty, sz, tz;

    theta0     =   CalculateTheta0(particle, *cross_sections.at(0)->GetMedium(), dr);

    rnd1 = SQRT2*theta0*erfInv( 2.*(RandomGenerator::Get().RandomDouble()-0.5) );
    rnd2 = SQRT2*theta0*erfInv( 2.*(RandomGenerator::Get().RandomDouble()-0.5) );

    sx      =   (rnd1/SQRT3+rnd2)/2;
    tx      =   rnd2;

    rnd1 = SQRT2*theta0*erfInv(2*(RandomGenerator::Get().RandomDouble()-0.5));
    rnd2 = SQRT2*theta0*erfInv(2*(RandomGenerator::Get().RandomDouble()-0.5));

    sy      =   (rnd1/SQRT3+rnd2)/2;
    ty      =   rnd2;

    sz      =   sqrt(max(1.-(sx*sx+sy*sy), 0.));
    tz      =   sqrt(max(1.-(tx*tx+ty*ty), 0.));


    Vector3D position;
    Vector3D direction;

    long double sinth, costh,sinph,cosph;
    sinth = (long double) sin(particle.GetDirection().GetTheta());
    costh = (long double) cos(particle.GetDirection().GetTheta());
    sinph = (long double) sin(particle.GetDirection().GetPhi());
    cosph = (long double) cos(particle.GetDirection().GetPhi());

    position = particle.GetPosition();

    // Rotation towards all tree axes
    direction = sz*particle.GetDirection();
    direction = direction + sx*Vector3D(costh*cosph, costh*sinph, -sinth);
    direction = direction + sy*Vector3D(-sinph, cosph, 0.);

    position = position + dr*direction;

    // Rotation towards all tree axes
    direction = tz*particle.GetDirection();
    direction = direction + tx*Vector3D(costh*cosph, costh*sinph, -sinth);
    direction = direction + ty*Vector3D(-sinph, cosph, 0.);

    direction.CalculateSphericalCoordinates();

    particle.SetPosition(position);
    particle.SetDirection(direction);
}

void ScatteringFirstOrder::EnableInterpolation(const PROPOSALParticle& particle,
                                            const std::vector<CrossSections*>& cross_sections,
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

