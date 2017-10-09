
#include <boost/bind.hpp>
#include <cmath>

#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/crossection/CrossSection.h"

#include "PROPOSAL/Constants.h"

#define UTILITY_INTEGRAL_IMPL(cls)                                                                                     \
    UtilityIntegral##cls::UtilityIntegral##cls(const Utility& utility)                                                 \
        : UtilityIntegral(utility)                                                                                     \
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

UtilityIntegral::UtilityIntegral(const UtilityIntegral& collection)
    : UtilityDecorator(collection)
    , integral_(collection.integral_)
{
}

UtilityIntegral::~UtilityIntegral()
{
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
        ei, ef, boost::bind(&UtilityIntegralDisplacement::FunctionToIntegral, this,  _1), 4, -rnd);
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

    for(std::vector<CrossSection*>::const_iterator iter = crosssections.begin(); iter != crosssections.end(); ++iter)
    {
        //TODO(mario): name Wed 2017/09/06
        // log_debug("Rate for %s = %f", crosssections_.at(i)->GetName().c_str(), rate);

        total_rate += (*iter)->CalculatedNdx(energy);
    }

    return UtilityDecorator::FunctionToIntegral(energy) * total_rate;
}

// ------------------------------------------------------------------------- //
double UtilityIntegralInteraction::Calculate(double ei, double ef, double rnd)
{
    (void) ef;

    return integral_.IntegrateWithRandomRatio(
        ei,
        ef,
        boost::bind(&UtilityIntegralInteraction::FunctionToIntegral, this,  _1),
        4,
        -rnd);
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

    //TODO(mario): Better way? Sat 2017/09/02
    double square_momentum = energy * energy - particle_def.mass * particle_def.mass;
    double particle_momentum = sqrt(std::max(square_momentum, 0.0));

    // return multiplier / max((particle_momentum / particle_.GetMass()) * particle_.GetLifetime() * SPEED, PARTICLE_POSITION_RESOLUTION);
    aux =  1.0 / std::max((particle_momentum / particle_def.mass) * particle_def.lifetime * SPEED, PARTICLE_POSITION_RESOLUTION);

    return UtilityDecorator::FunctionToIntegral(energy) * aux;
}

// ------------------------------------------------------------------------- //
double UtilityIntegralDecay::Calculate(double ei, double ef, double rnd)
{
    (void) ef;

    return integral_.IntegrateWithRandomRatio(
        ei,
        ef,
        boost::bind(&UtilityIntegralDecay::FunctionToIntegral, this,  _1),
        4,
        -rnd);
}

/******************************************************************************
*                            Utility Time                            *
******************************************************************************/

UTILITY_INTEGRAL_IMPL(Time)

// ------------------------------------------------------------------------- //
double UtilityIntegralTime::FunctionToIntegral(double energy)
{
    const ParticleDef& particle_def = utility_.GetParticleDef();
    //TODO(mario): Better way? Sat 2017/09/02
    double square_momentum = energy * energy - particle_def.mass * particle_def.mass;
    double particle_momentum = sqrt(std::max(square_momentum, 0.0));

    return energy  / (particle_momentum * SPEED) * UtilityDecorator::FunctionToIntegral(energy);
}

// ------------------------------------------------------------------------- //
double UtilityIntegralTime::Calculate(double ei, double ef, double rnd)
{
    (void) rnd;

    return integral_.Integrate(
        ei, ef, boost::bind(&UtilityIntegralTime::FunctionToIntegral, this,  _1), 4);
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

    for(std::vector<CrossSection*>::const_iterator iter = crosssections.begin(); iter != crosssections.end(); ++iter)
    {
        sum += (*iter)->CalculatedE2dx(energy);
    }

    return UtilityDecorator::FunctionToIntegral(energy) * sum;
}

// ------------------------------------------------------------------------- //
double UtilityIntegralContRand::Calculate(double ei, double ef, double rnd)
{
    (void) rnd;
    return integral_.Integrate(ei, ef, boost::bind(&UtilityIntegralContRand::FunctionToIntegral, this, _1), 4);
}

/******************************************************************************
*                             Utility Scattering                             *
******************************************************************************/

UTILITY_INTEGRAL_IMPL(Scattering)

// ------------------------------------------------------------------------- //
double UtilityIntegralScattering::FunctionToIntegral(double energy)
{
    double aux2;

    //TODO(mario): Better way? Sat 2017/09/02
    double square_momentum = energy * energy - utility_.GetParticleDef().mass * utility_.GetParticleDef().mass;
    aux2    =   RY* energy / square_momentum;

    return UtilityDecorator::FunctionToIntegral(energy) * aux2 * aux2;
}

// ------------------------------------------------------------------------- //
double UtilityIntegralScattering::Calculate(double ei, double ef, double rnd)
{
    (void) rnd;
    return integral_.Integrate(ei, ef, boost::bind(&UtilityIntegralScattering::FunctionToIntegral, this, _1), 4);
}

#undef UTILITY_INTEGRAL_DEC
