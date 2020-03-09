
#include <cmath>
#include <functional>

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"

#define UTILITY_INTEGRAL_IMPL(cls)                                             \
    UtilityIntegral##cls::UtilityIntegral##cls(CrossSectionList cross)         \
        : UtilityIntegral(cross)                                               \
        , displacement_(new UtilityIntegralDisplacement(cross))                \
    {                                                                          \
    }

using namespace PROPOSAL;

/******************************************************************************
 *                              Utility Integral                              *
 ******************************************************************************/

UtilityIntegral::UtilityIntegral(CrossSectionList cross)
    : UtilityDecorator(cross)
    , integral_(IROMB, IMAXS, IPREC2)
{
}

double UtilityIntegral::Calculate(double ei, double ef, double rnd)
{
    return integral_.IntegrateWithRandomRatio(ei, ef,
        std::bind(&UtilityIntegral::FunctionToIntegral, this, std::placeholders::_1), 4,
        -rnd);
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

    UtilityIntegralDisplacement::UtilityIntegralDisplacement(CrossSectionList cross)
        : UtilityIntegral(cross)
    {
    }

double UtilityIntegralDisplacement::FunctionToIntegral(double energy)
{
    double result = 0.0;
    for (const auto& cross : crosssections)
        result += cross->CalculatedEdx(energy);

    return -1.0 / result;
}

/******************************************************************************
 *                            Utility Interaction                            *
 ******************************************************************************/

UTILITY_INTEGRAL_IMPL(Interaction)
double UtilityIntegralInteraction::FunctionToIntegral(double energy)
{
    double total_rate = 0.0;
    for (const auto& cross : crosssections)
        total_rate += cross->CalculatedNdx(energy);

    return displacement_->FunctionToIntegral(energy) * total_rate;
}

/******************************************************************************
 *                            Utility Decay                            *
 ******************************************************************************/

UTILITY_INTEGRAL_IMPL(Decay)
double UtilityIntegralDecay::FunctionToIntegral(double energy)
{
    double lifetime = crosssections.front().GetParticleDef().lifetime;
    if (lifetime < 0) {
        return 0;
    }

    double mass = crosssections.front()->GetParametrization().GetParticleMass();
    double square_momentum = (energy - mass) * (energy + mass);
    double particle_momentum = std::sqrt(std::max(square_momentum, 0.0));

    double aux = 1.0
        / std::max((particle_momentum / mass) * lifetime * SPEED,
              PARTICLE_POSITION_RESOLUTION);

    return displacement_->FunctionToIntegral(energy) * aux;
}

/******************************************************************************
 *                            Utility Time                            *
 ******************************************************************************/

UTILITY_INTEGRAL_IMPL(Time)
double UtilityIntegralTime::FunctionToIntegral(double energy)
{
    double mass{ crosssections.front()->GetParametrization().GetParticleMass() };
    double square_momentum{ (energy - mass) * (energy + mass) };
    double particle_momentum{ std::sqrt(std::max(square_momentum, 0.0)) };

    return energy / (particle_momentum * SPEED)
        * displacement_->FunctionToIntegral(energy);
}

/******************************************************************************
 *                            Utility ContRand                            *
 ******************************************************************************/

UTILITY_INTEGRAL_IMPL(ContRand)
double UtilityIntegralContRand::FunctionToIntegral(double energy)
{
    double sum = 0.0;
    for (const auto& cross : crosssections)
        sum += cross->CalculatedE2dx(energy);

    return displacement_->FunctionToIntegral(energy) * sum;
}

/******************************************************************************
 *                             Utility Scattering                             *
 ******************************************************************************/

UTILITY_INTEGRAL_IMPL(Scattering)

// ------------------------------------------------------------------------- //
double UtilityIntegralScattering::FunctionToIntegral(double energy)
{
    double mass{ crosssections.front()->GetParametrization().GetParticleMass() };
    double square_momentum = (energy - mass) * (energy + mass);
    double aux = energy / square_momentum;

    return displacement_->FunctionToIntegral(energy) * aux * aux;
}

#undef UTILITY_INTEGRAL_DEC
