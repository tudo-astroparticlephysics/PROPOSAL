/*! \file   ScatteringHighland.cxx
 *   \brief  Source file for the ScatteringHighland routines.
 *
 *   For more details see the class documentation.
 *
 *   \date   2013.06.13
 *   \author Tomasz Fuchs
 */

#include <cmath>
// #include <algorithm>
// #include <stdlib.h>

#include "PROPOSAL/scattering/ScatteringHighland.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/particle/ParticleDef.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

ScatteringHighland::ScatteringHighland(const ParticleDef& particle_def, std::shared_ptr<const Medium> medium)
    : Scattering(particle_def), medium_(medium) {}

ScatteringHighland::ScatteringHighland(const ScatteringHighland& scattering)
    : Scattering(scattering), medium_(scattering.medium_) {}

ScatteringHighland::ScatteringHighland(const ParticleDef& particle_def,
                                       const ScatteringHighland& scattering)
    : Scattering(particle_def), medium_(scattering.medium_) {}

ScatteringHighland::~ScatteringHighland() {
}

bool ScatteringHighland::compare(const Scattering& scattering) const {
    const ScatteringHighland* scatteringHighland =
        dynamic_cast<const ScatteringHighland*>(&scattering);

    if (!scatteringHighland)
        return false;
    else if (*medium_ != *scatteringHighland->medium_)
        return false;
    else
        return true;
}

void ScatteringHighland::print(std::ostream& os) const
{
    os << "Medium:\n" << *medium_ << '\n';
}

// ------------------------------------------------------------------------- //
// Private methods
// ------------------------------------------------------------------------- //

double ScatteringHighland::CalculateTheta0(double dr, double ei, const Vector3D& pos) {
    // eq 6 of Lynch, Dahl
    // Nuclear Instruments and Methods in Physics Research Section B 58 (1991)
    double y = dr / medium_->GetRadiationLength(pos);
    double momentum_Sq = (ei - particle_def_.mass) * (ei + particle_def_.mass);
    double beta_p = momentum_Sq / ei; // beta * p = p^2/sqrt(p^2 + m^2)
    y = 13.6 * std::abs(particle_def_.charge) /
        (beta_p) * std::sqrt(y) * (1. + 0.088 * std::log10(y));
    return y;
}

//----------------------------------------------------------------------------//


Scattering::RandomAngles ScatteringHighland::CalculateRandomAngle(double dr,
                                                                  double ei,
                                                                  double ef,
                                                                  const Vector3D& pos,
                                                                  double rnd1,
                                                                  double rnd2,
                                                                  double rnd3,
                                                                  double rnd4) {
    (void)ef;

    double Theta0;
    Scattering::RandomAngles random_angles;

    Theta0 = CalculateTheta0(dr, ei, pos);

    rnd1 = Theta0 * inverseErrorFunction(rnd1);
    rnd2 = Theta0 * inverseErrorFunction(rnd2);

    random_angles.sx = 0.5 * (rnd1 / SQRT3 + rnd2);
    random_angles.tx = rnd2;

    rnd1 = Theta0 * inverseErrorFunction(rnd3);
    rnd2 = Theta0 * inverseErrorFunction(rnd4);

    random_angles.sy = 0.5 * (rnd1 / SQRT3 + rnd2);
    random_angles.ty = rnd2;

    return random_angles;
}
