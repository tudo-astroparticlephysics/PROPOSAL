
#include <boost/bind.hpp>

#include "PROPOSAL/math/MathModel.h"
#include "PROPOSAL/scattering/ScatteringDefault.h"
#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
// #include "PROPOSAL/crossection/CrossSections.h"
// #include "PROPOSAL/crossection/Ionization.h"
// #include "PROPOSAL/crossection/Bremsstrahlung.h"
// #include "PROPOSAL/crossection/Epairproduction.h"
// #include "PROPOSAL/crossection/Photonuclear.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;
using namespace std;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

ScatteringDefault::ScatteringDefault()
    : Scattering()
{
}

ScatteringDefault::ScatteringDefault(const ScatteringDefault &scattering)
    : Scattering(scattering)
{
}

//----------------------------------------------------------------------------//
long double ScatteringDefault::CalculateTheta0(const PROPOSALParticle& particle, const Medium& medium, double dr, double disp)
{

    double aux    = disp;
    double cutoff = 1;

    // TODO: Decide if using RadiationLength calculated in Bremsstrahlung or Medium
    double radiation_lenght = medium.GetRadiationLength();

    // TODO: check if one has to take the absolute value of the particle charge
    aux = sqrt(max(aux, 0.0) / radiation_lenght) * particle.GetCharge();
    aux *= max(1 + 0.038 * log(dr / radiation_lenght), 0.0);

    return min(aux, cutoff);
}

// ------------------------------------------------------------------------- //
Scattering::RandomAngles ScatteringDefault::CalculateRandomAngle(const PROPOSALParticle& particle, const Medium& medium, double dr, double disp)
{
    double Theta0,rnd1,rnd2;
    Scattering::RandomAngles random_angles;

    Theta0 = CalculateTheta0(particle, medium, dr, disp);

    rnd1 = SQRT2*Theta0*erfInv( 2.*(RandomGenerator::Get().RandomDouble()-0.5) );
    rnd2 = SQRT2*Theta0*erfInv( 2.*(RandomGenerator::Get().RandomDouble()-0.5) );

    random_angles.sx = (rnd1/SQRT3+rnd2)/2;
    random_angles.tx = rnd2;

    rnd1 = SQRT2*Theta0*erfInv(2*(RandomGenerator::Get().RandomDouble()-0.5));
    rnd2 = SQRT2*Theta0*erfInv(2*(RandomGenerator::Get().RandomDouble()-0.5));

    random_angles.sy = (rnd1/SQRT3+rnd2)/2;
    random_angles.ty = rnd2;

    return random_angles;
}
