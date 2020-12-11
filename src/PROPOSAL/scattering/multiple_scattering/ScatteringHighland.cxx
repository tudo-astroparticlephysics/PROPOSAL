/*! \file   ScatteringHighland.cxx
 *   \brief  Source file for the ScatteringHighland routines.
 *
 *   For more details see the class documentation.
 *
 *   \date   2013.06.13
 *   \author Tomasz Fuchs
 */

#include <cmath>
#include <array>
// #include <algorithm>
// #include <stdlib.h>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/scattering/multiple_scattering/ScatteringHighland.h"

using std::array;
using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

ScatteringHighland::ScatteringHighland(
    const ParticleDef& particle_def, Medium const& medium)
    : Scattering(particle_def)
    , medium_(medium)
    , charge(particle_def.charge)
{
}

ScatteringHighland::ScatteringHighland(const ScatteringHighland& scattering)
    : Scattering(scattering)
    , medium_(scattering.medium_)
    , charge(scattering.charge)
{
}

ScatteringHighland::ScatteringHighland(
    const ParticleDef& particle_def, const ScatteringHighland& scattering)
    : Scattering(particle_def)
    , medium_(scattering.medium_)
    , charge(scattering.charge)
{
}

bool ScatteringHighland::compare(const Scattering& scattering) const
{
    const ScatteringHighland* scatteringHighland
        = dynamic_cast<const ScatteringHighland*>(&scattering);

    if (!scatteringHighland)
        return false;
    else if (medium_ == scatteringHighland->medium_)
        return true;
    else
        return false;
}

void ScatteringHighland::print(std::ostream& os) const
{
    os << "Medium:\n" << medium_ << '\n';
}

// ------------------------------------------------------------------------- //
// Private methods
// ------------------------------------------------------------------------- //

double ScatteringHighland::CalculateTheta0(double grammage, double ei, double ef)
{
    (void) ef;

    // eq 6 of Lynch, Dahl
    // Nuclear Instruments and Methods in Physics Research Section B 58 (1991)
    double y = grammage / medium_.GetRadiationLength();
    double momentum_Sq = (ei - mass) * (ei + mass);
    double beta_p = momentum_Sq / ei; // beta * p = p^2/sqrt(p^2 + m^2)
    y = 13.6 * std::abs(charge) / (beta_p)*std::sqrt(y)
        * (1. + 0.088 * std::log10(y));
    return y;
}

//----------------------------------------------------------------------------//

Scattering::RandomAngles ScatteringHighland::CalculateRandomAngle(double grammage,
    double ei, double ef, const array<double, 4>& rnd)
{
    double Theta0;
    Scattering::RandomAngles random_angles;

    Theta0 = CalculateTheta0(grammage, ei, ef);

    auto rnd1 = Theta0 * normalppf(rnd[0]);
    auto rnd2 = Theta0 * normalppf(rnd[1]);

    random_angles.sx = 0.5 * (rnd1 / SQRT3 + rnd2);
    random_angles.tx = rnd2;

    rnd1 = Theta0 * normalppf(rnd[2]);
    rnd2 = Theta0 * normalppf(rnd[3]);

    random_angles.sy = 0.5 * (rnd1 / SQRT3 + rnd2);
    random_angles.ty = rnd2;

    return random_angles;
}
