
#include <array>
#include <cmath>

#include "PROPOSAL/Logging.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/scattering/ScatteringHighlandIntegral.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"

#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"

using std::array;
using namespace PROPOSAL;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

ScatteringHighlandIntegral::ScatteringHighlandIntegral(const ParticleDef& p_def,
    std::shared_ptr<const Medium> medium, CrossSectionList cross)
    : Scattering(p_def)
    , medium(medium)
    , scatter_(new UtilityIntegralScattering(cross, p_def))
{
}

ScatteringHighlandIntegral::ScatteringHighlandIntegral(const ParticleDef& p_def,
    std::shared_ptr<const Medium> medium,
    std::shared_ptr<InterpolationDef> interpolation_def, CrossSectionList cross)
    : Scattering(p_def)
    , medium(medium)
    , scatter_(
          new UtilityInterpolantScattering(cross, p_def, interpolation_def))
{
}

ScatteringHighlandIntegral::~ScatteringHighlandIntegral() {}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

bool ScatteringHighlandIntegral::compare(const Scattering& scattering) const
{
    const ScatteringHighlandIntegral* scatteringHighlandIntegral
        = dynamic_cast<const ScatteringHighlandIntegral*>(&scattering);

    if (!scatteringHighlandIntegral)
        return false;
    /* else if (*scatter_ != *scatteringHighlandIntegral->scatter_) */
    /*     return false; */
    else
        return true;
}

void ScatteringHighlandIntegral::print(std::ostream& os) const {}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//------------------------private member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
long double ScatteringHighlandIntegral::CalculateTheta0(
    double dr, double ei, double ef, const Vector3D& pos)
{
    double aux = scatter_->Calculate(ei, ef, 0.0)
        * medium->GetDensityDistribution().Evaluate(pos);
    double cutoff = 1;
    double radiation_length = medium->GetRadiationLength(pos);

    aux = 13.6 * std::sqrt(std::max(aux, 0.0) / radiation_length)
        * std::abs(particle_def_.charge);
    aux *= std::max(1 + 0.038 * std::log(dr / radiation_length), 0.0);

    return std::min(aux, cutoff);
}

//----------------------------------------------------------------------------//

Scattering::RandomAngles ScatteringHighlandIntegral::CalculateRandomAngle(
    double dr, double ei, double ef, const Vector3D& pos,
    const array<double, 4>& rnd)
{
    double Theta0;
    Scattering::RandomAngles random_angles;

    Theta0 = CalculateTheta0(dr, ei, ef, pos);

    auto rnd1 = Theta0 * inverseErrorFunction(rnd[0]);
    auto rnd2 = Theta0 * inverseErrorFunction(rnd[1]);

    random_angles.sx = 0.5 * (rnd1 / SQRT3 + rnd2);
    random_angles.tx = rnd2;

    rnd1 = Theta0 * inverseErrorFunction(rnd[2]);
    rnd2 = Theta0 * inverseErrorFunction(rnd[3]);

    random_angles.sy = 0.5 * (rnd1 / SQRT3 + rnd2);
    random_angles.ty = rnd2;

    return random_angles;
}
