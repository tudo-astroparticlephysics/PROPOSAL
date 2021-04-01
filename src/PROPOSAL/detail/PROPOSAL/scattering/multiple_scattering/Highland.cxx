/*! \file   ScatteringHighland.cxx
 *   \brief  Source file for the ScatteringHighland routines.
 *
 *   For more details see the class documentation.
 *
 *   \date   2013.06.13
 *   \author Tomasz Fuchs
 */

#include <array>
#include <cmath>
// #include <algorithm>
// #include <stdlib.h>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/scattering/multiple_scattering/Highland.h"

using namespace PROPOSAL::multiple_scattering;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

Highland::Highland(const ParticleDef& particle_def, Medium const& medium)
    : Parametrization(particle_def.mass)
    , charge(particle_def.charge)
    , radiation_length(medium.GetRadiationLength())
{
}

bool Highland::compare(const Parametrization& scattering) const
{
    auto sc = dynamic_cast<const Highland*>(&scattering);

    if (radiation_length != sc->radiation_length)
        return false;
    else if (charge != sc->charge)
        return false;
    return true;
}

void Highland::print(std::ostream& os) const
{
    os << "charge: " << charge << '\n'
       << "radiation_length: " << radiation_length << "\n";
}

// ------------------------------------------------------------------------- //
// Private methods
// ------------------------------------------------------------------------- //

double Highland::CalculateTheta0(double grammage, double ei, double ef)
{
    (void)ef;

    // eq 6 of Lynch, Dahl
    // Nuclear Instruments and Methods in Physics Research Section B 58 (1991)
    double y = grammage / radiation_length;
    double momentum_Sq = (ei - mass) * (ei + mass);
    double beta_p = momentum_Sq / ei; // beta * p = p^2/sqrt(p^2 + m^2)
    y = 13.6 * std::abs(charge) / (beta_p)*std::sqrt(y)
        * (1. + 0.088 * std::log10(y));
    return y;
}

//----------------------------------------------------------------------------//

ScatteringOffset Highland::CalculateRandomAngle(
    double grammage, double ei, double ef, const std::array<double, 4>& rnd)
{
    ScatteringOffset offsets;
    auto Theta0 = CalculateTheta0(grammage, ei, ef);

    auto rnd1 = Theta0 * normalppf(rnd[0]);
    auto rnd2 = Theta0 * normalppf(rnd[1]);

    offsets.sx = 0.5 * (rnd1 / SQRT3 + rnd2);
    offsets.tx = rnd2;

    rnd1 = Theta0 * normalppf(rnd[2]);
    rnd2 = Theta0 * normalppf(rnd[3]);

    offsets.sy = 0.5 * (rnd1 / SQRT3 + rnd2);
    offsets.ty = rnd2;

    return offsets;
}
