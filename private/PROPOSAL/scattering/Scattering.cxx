/*! \file   Scattering.cxx
*   \brief  Source filefor the Scattering bug routines.
*
*   This version has a major bug and produces too small scattering angles.
*
*   \date   2013.08.19
*   \author Tomasz Fuchs
**/

#include <cmath>

// #include "PROPOSAL/Output.h"
#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/particle/PROPOSALParticle.h"
#include "PROPOSAL/scattering/Scattering.h"

using namespace PROPOSAL;

/******************************************************************************
*                                 Scattering                                  *
******************************************************************************/


Scattering::Scattering(PROPOSALParticle& particle)
    : particle_(particle)
{
}

Scattering::Scattering(const Scattering& scattering)
    : particle_(scattering.particle_)
{
}

Scattering::~Scattering()
{
}

void Scattering::Scatter(double dr, double ei, double ef)
{
    double sz,tz;

    RandomAngles random_angles = CalculateRandomAngle(dr, ei, ef);

    sz = sqrt(std::max(1.-(random_angles.sx*random_angles.sx+random_angles.sy*random_angles.sy), 0.));
    tz = sqrt(std::max(1.-(random_angles.tx*random_angles.tx+random_angles.ty*random_angles.ty), 0.));

    Vector3D position;
    Vector3D direction;

    long double sinth, costh,sinph,cosph;
    sinth = (long double) sin(particle_.GetDirection().GetTheta());
    costh = (long double) cos(particle_.GetDirection().GetTheta());
    sinph = (long double) sin(particle_.GetDirection().GetPhi());
    cosph = (long double) cos(particle_.GetDirection().GetPhi());

    position = particle_.GetPosition();

    // Rotation towards all tree axes
    direction = sz*particle_.GetDirection();
    direction = direction + random_angles.sx*Vector3D(costh*cosph, costh*sinph, -sinth);
    direction = direction + random_angles.sy*Vector3D(-sinph, cosph, 0.);

    position = position + dr*direction;

    // Rotation towards all tree axes
    direction = tz*particle_.GetDirection();
    direction = direction + random_angles.tx*Vector3D(costh*cosph, costh*sinph, -sinth);
    direction = direction + random_angles.ty*Vector3D(-sinph, cosph, 0.);

    direction.CalculateSphericalCoordinates();

    particle_.SetPosition(position);
    particle_.SetDirection(direction);
}
