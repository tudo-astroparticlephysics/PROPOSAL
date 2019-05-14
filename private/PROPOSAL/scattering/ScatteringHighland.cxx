/*! \file   ScatteringHighland.cxx
 *   \brief  Source file for the ScatteringHighland routines.
 *
 *   For more details see the class documentation.
 *
 *   \date   2013.06.13
 *   \author Tomasz Fuchs
 */

// #include <cmath>
// #include <algorithm>
// #include <stdlib.h>

#include "PROPOSAL/scattering/ScatteringHighland.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/particle/Particle.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

ScatteringHighland::ScatteringHighland(Particle& particle, const Medium& medium)
    : Scattering(particle)
    , medium_(medium.clone())
{
}

ScatteringHighland::ScatteringHighland(const ScatteringHighland& scattering)
    : Scattering(scattering)
    , medium_(scattering.medium_->clone())
{
}

ScatteringHighland::ScatteringHighland(Particle& particle, const ScatteringHighland& scattering)
    : Scattering(particle)
    , medium_(scattering.medium_->clone())
{
}

ScatteringHighland::~ScatteringHighland()
{
    delete medium_;
}

bool ScatteringHighland::compare(const Scattering& scattering) const
{
    const ScatteringHighland* scatteringHighland = dynamic_cast<const ScatteringHighland*>(&scattering);

    if (!scatteringHighland)
        return false;
    else if (*medium_ != *scatteringHighland->medium_)
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
// Private methods
// ------------------------------------------------------------------------- //

double ScatteringHighland::CalculateTheta0(double dr)
{
    // eq 6 of Lynch, Dahl
    // Nuclear Instruments and Methods in Physics Research Section B 58 (1991)
    double y = dr / medium_->GetRadiationLength();
    double beta = 1. / sqrt(1 + particle_.GetMass() * particle_.GetMass() / (particle_.GetMomentum() * particle_.GetMomentum()));
    y = 13.6 * std::abs(particle_.GetCharge()) / (particle_.GetMomentum() * beta) * sqrt(y) * (1. + 0.088 * log10(y));
    return y;
}

//----------------------------------------------------------------------------//

Scattering::RandomAngles ScatteringHighland::CalculateRandomAngle(double dr, double ei, double ef)
{
    (void)ei;
    (void)ef;

    double Theta0, rnd1, rnd2;
    Scattering::RandomAngles random_angles;

    Theta0 = CalculateTheta0(dr);

    rnd1 = Theta0 * inverseErrorFunction(RandomGenerator::Get().RandomDouble());
    rnd2 = Theta0 * inverseErrorFunction(RandomGenerator::Get().RandomDouble());

    random_angles.sx = 0.5 * (rnd1 / SQRT3 + rnd2);
    random_angles.tx = rnd2;

    rnd1 = Theta0 * inverseErrorFunction(RandomGenerator::Get().RandomDouble());
    rnd2 = Theta0 * inverseErrorFunction(RandomGenerator::Get().RandomDouble());

    random_angles.sy = 0.5 * (rnd1 / SQRT3 + rnd2);
    random_angles.ty = rnd2;

    return random_angles;
}
