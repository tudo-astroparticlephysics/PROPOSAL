
#include <functional>

#include <cmath>

#include "PROPOSAL/crossection/IonizIntegral.h"
#include "PROPOSAL/crossection/parametrization/Ionization.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

IonizIntegral::IonizIntegral(const Ionization& param)
    : CrossSectionIntegral(DynamicData::DeltaE, param)
{
}

IonizIntegral::IonizIntegral(const IonizIntegral& brems)
    : CrossSectionIntegral(brems)
{
}

IonizIntegral::~IonizIntegral() {}

// ----------------------------------------------------------------- //
// Public methods
// ----------------------------------------------------------------- //

double IonizIntegral::CalculatedEdxWithoutMultiplier(double energy)
{
    double result, aux;

    Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);
    ParticleDef particle_def = parametrization_->GetParticleDef();
    const Medium& medium     = parametrization_->GetMedium();

    // TODO(mario): Better way? Sat 2017/09/02

    // PDG eq. 33.10
    // with Spin 1/2 correction by Rossi
    double square_momentum   = (energy - particle_def.mass) * (energy + particle_def.mass);
    double particle_momentum = std::sqrt(std::max(square_momentum, 0.0));
    double beta              = particle_momentum / energy;
    double gamma             = energy / particle_def.mass;

    aux    = beta * gamma / (1.e-6 * medium.GetI());
    result = std::log(limits.vUp * (2 * ME * energy)) + 2 * std::log(aux);
    aux    = limits.vUp / (2 * (1 + 1 / gamma));
    result += aux * aux;
    aux = beta * beta;
    result -= aux * (1 + limits.vUp / limits.vMax) + Delta(beta, gamma);

    if (result > 0)
    {
        result *= IONK * particle_def.charge * particle_def.charge * medium.GetZA() / (2 * aux);
    } else
    {
        result = 0;
    }
    return  medium.GetMassDensity() * result + energy * dedx_integral_.Integrate(
                        limits.vMin,
                        limits.vUp,
                        std::bind(&Parametrization::FunctionToDEdxIntegral, parametrization_, energy, std::placeholders::_1),
                        4);
}

double IonizIntegral::CalculatedEdx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    return parametrization_->GetMultiplier() * IonizIntegral::CalculatedEdxWithoutMultiplier(energy);
}

// ------------------------------------------------------------------------- //
double IonizIntegral::CalculatedE2dx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    return parametrization_->GetMultiplier() * IonizIntegral::CalculatedE2dxWithoutMultiplier(energy);
}

double IonizIntegral::CalculatedE2dxWithoutMultiplier(double energy)
{
    Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);

    return de2dx_integral_.Integrate(
        limits.vMin,
        limits.vUp,
        std::bind(&Parametrization::FunctionToDE2dxIntegral, parametrization_, energy, std::placeholders::_1),
        2);
}

// ------------------------------------------------------------------------- //
double IonizIntegral::CalculatedNdx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);
    ;
    sum_of_rates_ =
        dndx_integral_[0].Integrate(limits.vUp,
                                    limits.vMax,
                                    std::bind(&Parametrization::FunctionToDNdxIntegral, parametrization_, energy, std::placeholders::_1),
                                    3,
                                    1);

    return parametrization_->GetMultiplier() * sum_of_rates_;
}

// ------------------------------------------------------------------------- //
double IonizIntegral::CalculatedNdx(double energy, double rnd)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    rnd_ = rnd;

    Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);
    ;
    sum_of_rates_ = dndx_integral_[0].IntegrateWithRandomRatio(
        limits.vUp,
        limits.vMax,
        std::bind(&Parametrization::FunctionToDNdxIntegral, parametrization_, energy, std::placeholders::_1),
        3,
        rnd,
        1);

    return parametrization_->GetMultiplier() * sum_of_rates_;
}

// ------------------------------------------------------------------------- //
double IonizIntegral::CalculateStochasticLoss(double energy, double rnd1)
{
    double rnd, rsum;
    const Medium& medium = parametrization_->GetMedium();

    rnd  = medium.GetSumCharge() * rnd1;
    rsum = 0;

    for (unsigned int i = 0; i < components_.size(); i++)
    {
        rsum += components_[i]->GetAtomInMolecule() * components_[i]->GetNucCharge();

        if (rsum > rnd)
        {
            return energy * dndx_integral_[0].GetUpperLimit();
        }
    }

    log_fatal("SumCharge of medium was not initialized correctly");

    return 0;
}

// ------------------------------------------------------------------------- //
double IonizIntegral::Delta(double beta, double gamma)
{
    const Medium& medium = parametrization_->GetMedium();
    double X;

    X = std::log(beta * gamma) / std::log(10);

    if (X < medium.GetX0())
    {
        return medium.GetD0() * std::pow(10, 2 * (X - medium.GetX0()));
    } else if (X < medium.GetX1())
    {
        return 2 * LOG10 * X + medium.GetC() + medium.GetA() * std::pow(medium.GetX1() - X, medium.GetM());
    } else
    {
        return 2 * LOG10 * X + medium.GetC();
    }
}
