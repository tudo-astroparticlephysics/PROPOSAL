
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

using namespace PROPOSAL;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

ScatteringHighlandIntegral::ScatteringHighlandIntegral(const ParticleDef& particle_def, const Utility& utility)
    : Scattering(particle_def)
    , scatter_(new UtilityIntegralScattering(utility))
{
    if (particle_def != utility.GetParticleDef())
    {
        log_fatal("Particle definition should be equal to the utility particle definition!");
    }
}

ScatteringHighlandIntegral::ScatteringHighlandIntegral(const ParticleDef& particle_def,
                                                       const Utility& utility,
                                                       const InterpolationDef& interpolation_def)
    : Scattering(particle_def)
    , scatter_(new UtilityInterpolantScattering(utility, interpolation_def))
{
    if (particle_def != utility.GetParticleDef())
    {
        log_fatal("Particle definition should be equal to the utility particle definition!");
    }
}

ScatteringHighlandIntegral::ScatteringHighlandIntegral(const ScatteringHighlandIntegral& scattering)
    : Scattering(scattering)
    , scatter_(scattering.scatter_->clone(scattering.scatter_->GetUtility()))
{
}

ScatteringHighlandIntegral::ScatteringHighlandIntegral(const ParticleDef& particle_def,
                                                       const Utility& utility,
                                                       const ScatteringHighlandIntegral& scattering)
    : Scattering(particle_def)
    , scatter_(scattering.scatter_->clone(utility))
{
    if (particle_def != scattering.GetParticleDef())
    {
        log_fatal("Particle definition should be equal to the scattering particle definition!");
    }
    if (utility != scattering.scatter_->GetUtility())
    {
        log_fatal("Utilities of the ScatteringHighlandIntegral should have same values!");
    }
}

ScatteringHighlandIntegral::~ScatteringHighlandIntegral()
{
    delete scatter_;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

bool ScatteringHighlandIntegral::compare(const Scattering& scattering) const
{
    const ScatteringHighlandIntegral* scatteringHighlandIntegral =
        dynamic_cast<const ScatteringHighlandIntegral*>(&scattering);

    if (!scatteringHighlandIntegral)
        return false;
    else if (*scatter_ != *scatteringHighlandIntegral->scatter_)
        return false;
    else
        return true;
}

void ScatteringHighlandIntegral::print(std::ostream& os) const
{
    os << "Propagation Utility:\n" << scatter_->GetUtility() << '\n';
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//------------------------private member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
long double ScatteringHighlandIntegral::CalculateTheta0(double dr, double ei, double ef, const Vector3D& pos)
{
    double aux              = scatter_->Calculate(ei, ef, 0.0) *
                              scatter_->GetUtility().GetMedium()->GetDensityDistribution().Evaluate(pos);
    double cutoff           = 1;
    double radiation_length = scatter_->GetUtility().GetMedium()->GetRadiationLength(pos);

    aux = 13.6 * std::sqrt(std::max(aux, 0.0) / radiation_length) * std::abs(particle_def_.charge);
    aux *= std::max(1 + 0.038 * std::log(dr / radiation_length), 0.0);

    return std::min(aux, cutoff);
}

//----------------------------------------------------------------------------//

Scattering::RandomAngles ScatteringHighlandIntegral::CalculateRandomAngle(double dr,
                                                                          double ei,
                                                                          double ef,
                                                                          const Vector3D& pos,
                                                                          double rnd1,
                                                                          double rnd2,
                                                                          double rnd3,
                                                                          double rnd4)
{
    double Theta0;
    Scattering::RandomAngles random_angles;

    Theta0 = CalculateTheta0(dr, ei, ef, pos);

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
