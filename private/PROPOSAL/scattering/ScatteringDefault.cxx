
#include <boost/bind.hpp>

#include "PROPOSAL/math/MathModel.h"
#include "PROPOSAL/scattering/ScatteringDefault.h"
#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Constants.h"

#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"


using namespace PROPOSAL;
using namespace std;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

ScatteringDefault::ScatteringDefault(PROPOSALParticle& particle, Utility& utility)
    : Scattering(particle)
    , scatter(new UtilityIntegralScattering(utility))
{
}

ScatteringDefault::ScatteringDefault(PROPOSALParticle& particle, Utility& utility, InterpolationDef interpolation_def)
    : Scattering(particle)
    , scatter(new UtilityInterpolantScattering(utility, interpolation_def))
{
}

ScatteringDefault::ScatteringDefault(const ScatteringDefault& scattering)
    : Scattering(scattering)
    , scatter(scattering.scatter->clone())
{
}

ScatteringDefault::ScatteringDefault(PROPOSALParticle& particle, const ScatteringDefault& scattering)
    : Scattering(particle)
    , scatter(scattering.scatter->clone())
{
}

ScatteringDefault::~ScatteringDefault()
{
    delete scatter;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// ScatteringDefault& ScatteringDefault::operator=(const ScatteringDefault &scattering){
//     if (this != &scattering)
//     {
//       ScatteringDefault tmp(scattering);
//       swap(tmp);
//     }
//     return *this;
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// bool ScatteringDefault::operator==(const ScatteringDefault &scattering) const
// {
//     if( interpolant_ != NULL && scattering.interpolant_ != NULL)
//     {
//         if( *interpolant_   != *scattering.interpolant_) return false;
//     }
//
//     if( interpolant_diff_ != NULL && scattering.interpolant_diff_ != NULL)
//     {
//         if( *interpolant_diff_   != *scattering.interpolant_diff_) return false;
//     }
//
//     if( integral_   != scattering.integral_) return false;
//
//     //else
//     return true;
//
// }
//
//
// //----------------------------------------------------------------------------//
// //----------------------------------------------------------------------------//
//
//
// bool ScatteringDefault::operator!=(const ScatteringDefault &scattering) const {
//   return !(*this == scattering);
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// void ScatteringDefault::swap(ScatteringDefault &scattering)
// {
//     using std::swap;
//
//     swap(do_interpolation_,scattering.do_interpolation_);
//
//     if(scattering.interpolant_ != NULL)
//     {
//         interpolant_->swap(*scattering.interpolant_);
//     }
//     else
//     {
//         interpolant_ = NULL;
//     }
//
//     if(scattering.interpolant_diff_ != NULL)
//     {
//         interpolant_diff_->swap(*scattering.interpolant_diff_) ;
//     }
//     else
//     {
//         interpolant_diff_ = NULL;
//     }
//
//     integral_.swap(scattering.integral_);
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//------------------------private member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
long double ScatteringDefault::CalculateTheta0(double dr, double ei, double ef)
{
    double aux              = scatter->Calculate(ei, ef, 0.0);
    double cutoff           = 1;
    double radiation_lenght = scatter->GetUtility().GetMedium().GetRadiationLength();

    // TODO: check if one has to take the absolute value of the particle charge
    aux = sqrt(max(aux, 0.0) / radiation_lenght) * particle_.GetCharge();
    aux *= max(1 + 0.038 * log(dr / radiation_lenght), 0.0);

    return min(aux, cutoff);
}

//----------------------------------------------------------------------------//

Scattering::RandomAngles ScatteringDefault::CalculateRandomAngle(double dr, double ei, double ef)
{
    double Theta0, rnd1, rnd2;
    Scattering::RandomAngles random_angles;

    Theta0 = CalculateTheta0(dr, ei, ef);

    rnd1 = SQRT2 * Theta0 * erfInv(2. * (RandomGenerator::Get().RandomDouble() - 0.5));
    rnd2 = SQRT2 * Theta0 * erfInv(2. * (RandomGenerator::Get().RandomDouble() - 0.5));

    random_angles.sx = (rnd1 / SQRT3 + rnd2) / 2;
    random_angles.tx = rnd2;

    rnd1 = SQRT2 * Theta0 * erfInv(2 * (RandomGenerator::Get().RandomDouble() - 0.5));
    rnd2 = SQRT2 * Theta0 * erfInv(2 * (RandomGenerator::Get().RandomDouble() - 0.5));

    random_angles.sy = (rnd1 / SQRT3 + rnd2) / 2;
    random_angles.ty = rnd2;

    return random_angles;
}
