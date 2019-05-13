
#include <functional>
#include <cmath>

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"

#define UTILITY_INTEGRAL_IMPL(cls)                                                                                     \
    UtilityIntegral##cls::UtilityIntegral##cls(const Utility& utility)                                                 \
        : UtilityIntegral(utility)                                                                                     \
    {                                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    UtilityIntegral##cls::UtilityIntegral##cls(const Utility& utility, const UtilityIntegral##cls& decorator)          \
        : UtilityIntegral(utility, decorator)                                                                          \
    {                                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    UtilityIntegral##cls::UtilityIntegral##cls(const UtilityIntegral##cls& decorator)                                  \
        : UtilityIntegral(decorator.utility_)                                                                          \
    {                                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    UtilityIntegral##cls::~UtilityIntegral##cls() {}

using namespace PROPOSAL;

/******************************************************************************
 *                              Utility Integral                              *
 ******************************************************************************/

UtilityIntegral::UtilityIntegral(const Utility& utility)
    : UtilityDecorator(utility)
    , integral_(IROMB, IMAXS, IPREC2)
{
}

UtilityIntegral::UtilityIntegral(const Utility& utility, const UtilityIntegral& collection)
    : UtilityDecorator(utility)
    , integral_(collection.integral_)
{
    if (utility != collection.GetUtility())
    {
        log_fatal("Utilities of the decorators should have same values!");
    }
}

UtilityIntegral::UtilityIntegral(const UtilityIntegral& collection)
    : UtilityDecorator(collection)
    , integral_(collection.integral_)
{
}

UtilityIntegral::~UtilityIntegral() {}

bool UtilityIntegral::compare(const UtilityDecorator& utility_decorator) const
{
    const UtilityIntegral* utility_integral = static_cast<const UtilityIntegral*>(&utility_decorator);

    if (integral_ != utility_integral->integral_)
        return false;
    else
        return true;
}

double UtilityIntegral::GetUpperLimit(double ei, double rnd)
{
    (void)ei;
    (void)rnd;

    return integral_.GetUpperLimit();
}

/******************************************************************************
 *                            Utility Displacement                            *
 ******************************************************************************/

UTILITY_INTEGRAL_IMPL(Displacement)

// ------------------------------------------------------------------------- //
double UtilityIntegralDisplacement::Calculate(double ei, double ef, double rnd)
{
    return integral_.IntegrateWithRandomRatio(
        ei, ef, std::bind(&UtilityIntegralDisplacement::FunctionToIntegral, this, std::placeholders::_1), 4, -rnd);
}

/******************************************************************************
 *                            Utility Interaction                            *
 ******************************************************************************/

UTILITY_INTEGRAL_IMPL(Interaction)

// ------------------------------------------------------------------------- //
double UtilityIntegralInteraction::FunctionToIntegral(double energy)
{
    double total_rate = 0.0;

    const std::vector<CrossSection*>& crosssections = utility_.GetCrosssections();

    for (std::vector<CrossSection*>::const_iterator iter = crosssections.begin(); iter != crosssections.end(); ++iter)
    {
        total_rate += (*iter)->CalculatedNdx(energy);
    }

    return UtilityDecorator::FunctionToIntegral(energy) * total_rate;
}

// ------------------------------------------------------------------------- //
double UtilityIntegralInteraction::Calculate(double ei, double ef, double rnd)
{
    (void)ef;

    return integral_.IntegrateWithRandomRatio(
        ei, ef, std::bind(&UtilityIntegralInteraction::FunctionToIntegral, this, std::placeholders::_1), 4, -rnd);
}

/******************************************************************************
 *                            Utility Decay                            *
 ******************************************************************************/

UTILITY_INTEGRAL_IMPL(Decay)

// ------------------------------------------------------------------------- //
double UtilityIntegralDecay::FunctionToIntegral(double energy)
{
    const ParticleDef& particle_def = utility_.GetParticleDef();
    double aux;

    if (particle_def.lifetime < 0)
    {
        return 0;
    }

    // TODO(mario): Better way? Sat 2017/09/02
    double square_momentum   = (energy - particle_def.mass) * (energy + particle_def.mass);
    double particle_momentum = std::sqrt(std::max(square_momentum, 0.0));

    aux = 1.0 / std::max((particle_momentum / particle_def.mass) * particle_def.lifetime * SPEED,
                         PARTICLE_POSITION_RESOLUTION);

    return UtilityDecorator::FunctionToIntegral(energy) * aux;
}

// ------------------------------------------------------------------------- //
double UtilityIntegralDecay::Calculate(double ei, double ef, double rnd)
{
    (void)ef;

    return integral_.IntegrateWithRandomRatio(
        ei, ef, std::bind(&UtilityIntegralDecay::FunctionToIntegral, this, std::placeholders::_1), 4, -rnd);
}

/******************************************************************************
 *                            Utility Time                            *
 ******************************************************************************/

UTILITY_INTEGRAL_IMPL(Time)

// ------------------------------------------------------------------------- //
double UtilityIntegralTime::FunctionToIntegral(double energy)
{
    // TODO(mario): Better way? Sat 2017/09/02
    double square_momentum   = (energy - utility_.GetParticleDef().mass) * (energy + utility_.GetParticleDef().mass);
    double particle_momentum = std::sqrt(std::max(square_momentum, 0.0));

    return energy / (particle_momentum * SPEED) * UtilityDecorator::FunctionToIntegral(energy);
}

// ------------------------------------------------------------------------- //
double UtilityIntegralTime::Calculate(double ei, double ef, double rnd)
{
    (void)rnd;

    return integral_.Integrate(ei, ef, std::bind(&UtilityIntegralTime::FunctionToIntegral, this, std::placeholders::_1), 4);
}

/******************************************************************************
 *                            Utility ContRand                            *
 ******************************************************************************/

UTILITY_INTEGRAL_IMPL(ContRand)

// ------------------------------------------------------------------------- //
double UtilityIntegralContRand::FunctionToIntegral(double energy)
{
    double sum = 0.0;

    const std::vector<CrossSection*>& crosssections = utility_.GetCrosssections();

    for (std::vector<CrossSection*>::const_iterator iter = crosssections.begin(); iter != crosssections.end(); ++iter)
    {
        sum += (*iter)->CalculatedE2dx(energy);
    }

    return UtilityDecorator::FunctionToIntegral(energy) * sum;
}

// ------------------------------------------------------------------------- //
double UtilityIntegralContRand::Calculate(double ei, double ef, double rnd)
{
    (void)rnd;
    return integral_.Integrate(ei, ef, std::bind(&UtilityIntegralContRand::FunctionToIntegral, this, std::placeholders::_1), 4);
}

/******************************************************************************
 *                             Utility Scattering                             *
 ******************************************************************************/

UTILITY_INTEGRAL_IMPL(Scattering)

// ------------------------------------------------------------------------- //
double UtilityIntegralScattering::FunctionToIntegral(double energy)
{
    double aux2;

    // TODO(mario): Better way? Sat 2017/09/02
    double square_momentum = (energy - utility_.GetParticleDef().mass) * (energy + utility_.GetParticleDef().mass);
    aux2                   = energy / square_momentum;

    return UtilityDecorator::FunctionToIntegral(energy) * aux2 * aux2;
}

// ------------------------------------------------------------------------- //
double UtilityIntegralScattering::Calculate(double ei, double ef, double rnd)
{
    (void)rnd;
    return integral_.Integrate(ei, ef, std::bind(&UtilityIntegralScattering::FunctionToIntegral, this, std::placeholders::_1), 4);
}

#undef UTILITY_INTEGRAL_DEC
