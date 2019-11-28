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
#include "PROPOSAL/math/RandomGenerator.h"

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

Directions Scattering::Scatter(double dr, double ei, double ef)
{
    double rnd1 = RandomGenerator::Get().RandomDouble();
    double rnd2 = RandomGenerator::Get().RandomDouble();
    double rnd3 = RandomGenerator::Get().RandomDouble();
    double rnd4 = RandomGenerator::Get().RandomDouble();

    return Scattering::Scatter(dr, ei, ef, rnd1, rnd2, rnd3, rnd4);
}

Directions Scattering::Scatter(double dr, double ei, double ef, double rnd1, double rnd2, double rnd3, double rnd4)
{
    Directions directions_;
    // u averaged continous propagation direction
    // n_i direction after continous propagation
    directions_.n_i_ = Vector3D(particle_.GetDirection());

    if(dr<=0)
    {
        return directions_;
    }

    double sz, tz;

    RandomAngles random_angles = CalculateRandomAngle(dr, ei, ef, rnd1, rnd2, rnd3, rnd4);

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
    directions_.u_ = sz * old_direction;
    directions_.u_ = directions_.u_ + random_angles.sx * rotate_vector_x;
    directions_.u_ = directions_.u_ + random_angles.sy * rotate_vector_y;

    // Rotation towards all tree axes
    directions_.n_i_ = tz * old_direction;
    directions_.n_i_ = directions_.n_i_ + random_angles.tx * rotate_vector_x;
    directions_.n_i_ = directions_.n_i_ + random_angles.ty * rotate_vector_y;
    directions_.n_i_.CalculateSphericalCoordinates();

    return directions_;
}

Scattering::RandomAngles Scattering::CalculateRandomAngle(double dr, double ei, double ef)
{
    double rnd1 = RandomGenerator::Get().RandomDouble();
    double rnd2 = RandomGenerator::Get().RandomDouble();
    double rnd3 = RandomGenerator::Get().RandomDouble();
    double rnd4 = RandomGenerator::Get().RandomDouble();

    return this->CalculateRandomAngle(dr, ei, ef, rnd1, rnd2, rnd3, rnd4);
}
