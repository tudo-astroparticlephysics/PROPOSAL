

#include "PROPOSAL/Logging.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/scattering/ScatteringHighlandIntegral.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/Particle.h"

#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"

using namespace PROPOSAL;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

ScatteringHighlandIntegral::ScatteringHighlandIntegral(Particle& particle, const Utility& utility)
    : Scattering(particle)
    , scatter_(new UtilityIntegralScattering(utility))
{
    if (particle.GetParticleDef() != utility.GetParticleDef())
    {
        log_fatal("Particle definition should be equal to the utility paricle definition!");
    }
}

ScatteringHighlandIntegral::ScatteringHighlandIntegral(Particle& particle,
                                                       const Utility& utility,
                                                       const InterpolationDef& interpolation_def)
    : Scattering(particle)
    , scatter_(new UtilityInterpolantScattering(utility, interpolation_def))
{
    if (particle.GetParticleDef() != utility.GetParticleDef())
    {
        log_fatal("Particle definition should be equal to the utility paricle definition!");
    }
}

ScatteringHighlandIntegral::ScatteringHighlandIntegral(const ScatteringHighlandIntegral& scattering)
    : Scattering(scattering)
    , scatter_(scattering.scatter_->clone(scattering.scatter_->GetUtility()))
{
}

ScatteringHighlandIntegral::ScatteringHighlandIntegral(Particle& particle,
                                                       const Utility& utility,
                                                       const ScatteringHighlandIntegral& scattering)
    : Scattering(particle)
    , scatter_(scattering.scatter_->clone(utility))
{
    if (particle.GetParticleDef() != scattering.GetParticle().GetParticleDef())
    {
        log_fatal("Particle definition should be equal to the scattering paricle definition!");
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

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//------------------------private member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
long double ScatteringHighlandIntegral::CalculateTheta0(double dr, double ei, double ef)
{
    double aux              = scatter_->Calculate(ei, ef, 0.0);
    double cutoff           = 1;
    double radiation_lenght = scatter_->GetUtility().GetMedium().GetRadiationLength();

    aux = 13.6 * std::sqrt(std::max(aux, 0.0) / radiation_lenght) * std::abs(particle_.GetCharge());
    aux *= std::max(1 + 0.038 * std::log(dr / radiation_lenght), 0.0);

    return std::min(aux, cutoff);
}

//----------------------------------------------------------------------------//

Scattering::RandomAngles ScatteringHighlandIntegral::CalculateRandomAngle(double dr, double ei, double ef)
{
    double Theta0, rnd1, rnd2;
    Scattering::RandomAngles random_angles;

    Theta0 = CalculateTheta0(dr, ei, ef);

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
