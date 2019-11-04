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

std::pair<Vector3D, Vector3D> Scattering::Scatter(double dr, double ei, double ef)
{
    std::pair<Vector3D, Vector3D> new_direction(Vector3D(0,0,0), particle_.GetDirection());

    if(dr<=0)
    {
        return new_direction;
    }

    double sz, tz;

    RandomAngles random_angles = CalculateRandomAngle(dr, ei, ef);

    sz = std::sqrt(std::max(1. - (random_angles.sx * random_angles.sx + random_angles.sy * random_angles.sy), 0.));
    tz = std::sqrt(std::max(1. - (random_angles.tx * random_angles.tx + random_angles.ty * random_angles.ty), 0.));

    const Vector3D old_direction = particle_.GetDirection();

    long double sinth, costh, sinph, cosph;
    sinth = (long double)std::sin(old_direction.GetTheta());
    costh = (long double)std::cos(old_direction.GetTheta());
    sinph = (long double)std::sin(old_direction.GetPhi());
    cosph = (long double)std::cos(old_direction.GetPhi());

    const Vector3D rotate_vector_x = Vector3D(costh * cosph, costh * sinph, -sinth);
    const Vector3D rotate_vector_y = Vector3D(-sinph, cosph, 0.);

    // Rotation towards all tree axes
    new_direction.first = sz * old_direction;
    new_direction.first = new_direction.first + random_angles.sx * rotate_vector_x;
    new_direction.first = new_direction.first + random_angles.sy * rotate_vector_y;

    // Rotation towards all tree axes
    new_direction.second = tz * old_direction;
    new_direction.second = new_direction.second + random_angles.tx * rotate_vector_x;
    new_direction.second = new_direction.second + random_angles.ty * rotate_vector_y;
    new_direction.second.CalculateSphericalCoordinates();

    return new_direction;
}
