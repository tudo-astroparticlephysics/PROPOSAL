/*! \file   Scattering.cxx
 *   \brief  Source filefor the Scattering bug routines.
 *
 *   This version has a major bug and produces too small scattering angles.
 *
 *   \date   2013.08.19
 *   \author Tomasz Fuchs
 **/

#include <cmath>

#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/scattering/Scattering.h"

using namespace PROPOSAL;

/******************************************************************************
 *                                 Scattering                                  *
 ******************************************************************************/

Scattering::Scattering(Particle& particle)
    : particle_(particle)
{
}

Scattering::Scattering(const Scattering& scattering)
    : particle_(scattering.particle_)
{
}

Scattering::~Scattering() {}

bool Scattering::operator==(const Scattering& scattering) const
{
    if (particle_ != scattering.particle_)
        return false;
    else
        return this->compare(scattering);
}

bool Scattering::operator!=(const Scattering& scattering) const
{
    return !(*this == scattering);
}

void Scattering::Scatter(double dr, double ei, double ef)
{
    double sz, tz;

    RandomAngles random_angles = CalculateRandomAngle(dr, ei, ef);

    sz = std::sqrt(std::max(1. - (random_angles.sx * random_angles.sx + random_angles.sy * random_angles.sy), 0.));
    tz = std::sqrt(std::max(1. - (random_angles.tx * random_angles.tx + random_angles.ty * random_angles.ty), 0.));

    Vector3D position;
    Vector3D direction;
    const Vector3D old_direction = particle_.GetDirection();

    long double sinth, costh, sinph, cosph;
    sinth = (long double)std::sin(old_direction.GetTheta());
    costh = (long double)std::cos(old_direction.GetTheta());
    sinph = (long double)std::sin(old_direction.GetPhi());
    cosph = (long double)std::cos(old_direction.GetPhi());

    const Vector3D rotate_vector_x = Vector3D(costh * cosph, costh * sinph, -sinth);
    const Vector3D rotate_vector_y = Vector3D(-sinph, cosph, 0.);

    position = particle_.GetPosition();

    // Rotation towards all tree axes
    direction = sz * old_direction;
    direction = direction + random_angles.sx * rotate_vector_x;
    direction = direction + random_angles.sy * rotate_vector_y;

    position = position + dr * direction;

    // Rotation towards all tree axes
    direction = tz * old_direction;
    direction = direction + random_angles.tx * rotate_vector_x;
    direction = direction + random_angles.ty * rotate_vector_y;

    direction.CalculateSphericalCoordinates();

    particle_.SetPosition(position);
    particle_.SetDirection(direction);
}
