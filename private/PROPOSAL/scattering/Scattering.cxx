/*! \file   Scattering.cxx
 *   \brief  Source filefor the Scattering bug routines.
 *
 *   This version has a major bug and produces too small scattering angles.
 *
 *   \date   2013.08.19
 *   \author Tomasz Fuchs
 **/

#include <cmath>
#include <memory>

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

std::shared_ptr<std::pair<Vector3D, Vector3D>> Scattering::Scatter(double dr, double ei, double ef)
{
    Vector3D u(0,0,0); // averaged continous propagation direction
    Vector3D n_i(particle_.GetDirection()); // direction after continous propagation

    if(dr<=0)
    {
        return std::make_shared<std::pair<Vector3D, Vector3D>>(std::make_pair(u, n_i));
    }

    double sz, tz;

    RandomAngles random_angles = CalculateRandomAngle(dr, ei, ef);

    sz = std::sqrt(std::max(1. - (random_angles.sx * random_angles.sx + random_angles.sy * random_angles.sy), 0.));
    tz = std::sqrt(std::max(1. - (random_angles.tx * random_angles.tx + random_angles.ty * random_angles.ty), 0.));

    const Vector3D old_direction = particle_.GetDirection();

    double sinth, costh, sinph, cosph;
    sinth = std::sin(old_direction.GetTheta());
    costh = std::cos(old_direction.GetTheta());
    sinph = std::sin(old_direction.GetPhi());
    cosph = std::cos(old_direction.GetPhi());

    const Vector3D rotate_vector_x = Vector3D(costh * cosph, costh * sinph, -sinth);
    const Vector3D rotate_vector_y = Vector3D(-sinph, cosph, 0.);

    // Rotation towards all tree axes
    u = sz * old_direction;
    u = u + random_angles.sx * rotate_vector_x;
    u = u + random_angles.sy * rotate_vector_y;

    // Rotation towards all tree axes
    n_i = tz * old_direction;
    n_i = n_i + random_angles.tx * rotate_vector_x;
    n_i = n_i + random_angles.ty * rotate_vector_y;
    n_i.CalculateSphericalCoordinates();

    return std::make_shared<std::pair<Vector3D, Vector3D>>(std::make_pair(u, n_i));
}
