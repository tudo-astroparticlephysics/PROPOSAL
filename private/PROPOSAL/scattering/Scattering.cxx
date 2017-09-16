/*! \file   Scattering.cxx
*   \brief  Source filefor the Scattering bug routines.
*
*   This version has a major bug and produces too small scattering angles.
*
*   \date   2013.08.19
*   \author Tomasz Fuchs
**/

#include <cmath>

#include "PROPOSAL/Output.h"
#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/particle/PROPOSALParticle.h"
#include "PROPOSAL/scattering/Scattering.h"
#include "PROPOSAL/scattering/ScatteringDefault.h"
#include "PROPOSAL/scattering/ScatteringMoliere.h"
#include "PROPOSAL/scattering/ScatteringFirstOrder.h"

using namespace PROPOSAL;

/******************************************************************************
*                                 Scattering                                  *
******************************************************************************/


Scattering::Scattering()
{
}

Scattering::~Scattering()
{
}

void Scattering::Scatter(PROPOSALParticle& particle, const Medium& medium, double dr, double disp)
{
    double sz,tz;

    RandomAngles random_angles = CalculateRandomAngle(particle, medium, dr, disp);

    sz = sqrt(std::max(1.-(random_angles.sx*random_angles.sx+random_angles.sy*random_angles.sy), 0.));
    tz = sqrt(std::max(1.-(random_angles.tx*random_angles.tx+random_angles.ty*random_angles.ty), 0.));

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
    direction = direction + random_angles.sx*Vector3D(costh*cosph, costh*sinph, -sinth);
    direction = direction + random_angles.sy*Vector3D(-sinph, cosph, 0.);

    position = position + dr*direction;

    // Rotation towards all tree axes
    direction = tz*particle.GetDirection();
    direction = direction + random_angles.tx*Vector3D(costh*cosph, costh*sinph, -sinth);
    direction = direction + random_angles.ty*Vector3D(-sinph, cosph, 0.);

    direction.CalculateSphericalCoordinates();

    particle.SetPosition(position);
    particle.SetDirection(direction);
}
