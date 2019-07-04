
#include <functional>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/EpairIntegral.h"
#include "PROPOSAL/crossection/parametrization/EpairProduction.h"
#include "PROPOSAL/medium/Medium.h"

using namespace PROPOSAL;

EpairIntegral::EpairIntegral(const EpairProduction& param)
    : CrossSectionIntegral(DynamicData::Epair, param)
{
}

EpairIntegral::EpairIntegral(const EpairIntegral& brems)
    : CrossSectionIntegral(brems)
{
}

EpairIntegral::~EpairIntegral() {}

// ----------------------------------------------------------------- //
// Public methods
// ----------------------------------------------------------------- //

double EpairIntegral::CalculatedEdx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    return parametrization_->GetMultiplier() * EpairIntegral::CalculatedEdxWithoutMultiplier(energy);
}

double EpairIntegral::CalculatedEdxWithoutMultiplier(double energy)
{
    double sum = 0;

    for (int i = 0; i < parametrization_->GetMedium().GetNumComponents(); i++)
    {
        parametrization_->SetCurrentComponent(i);
        Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);

        double r1  = 0.8;
        double rUp = limits.vUp * (1 - HALF_PRECISION);
        bool rflag = false;

        if (r1 < rUp)
        {
            if (2 * parametrization_->FunctionToDEdxIntegral(energy, r1) <
                parametrization_->FunctionToDEdxIntegral(energy, rUp))
            {
                rflag = true;
            }
        }

        if (rflag)
        {
            if (r1 > limits.vUp)
            {
                r1 = limits.vUp;
            }

            if (r1 < limits.vMin)
            {
                r1 = limits.vMin;
            }

            sum += dedx_integral_.Integrate(
                limits.vMin,
                r1,
                std::bind(&Parametrization::FunctionToDEdxIntegral, parametrization_, energy, std::placeholders::_1),
                4);
            double r2 = std::max(1 - limits.vUp, COMPUTER_PRECISION);

            if (r2 > 1 - r1)
            {
                r2 = 1 - r1;
            }

            sum +=
                dedx_integral_.Integrate(1 - limits.vUp,
                                         r2,
                                         std::bind(&EpairIntegral::FunctionToDEdxIntegralReverse, this, energy, std::placeholders::_1),
                                         2) +
                dedx_integral_.Integrate(
                    r2, 1 - r1, std::bind(&EpairIntegral::FunctionToDEdxIntegralReverse, this, energy, std::placeholders::_1), 4);

        }

        else
        {
            sum += dedx_integral_.Integrate(
                limits.vMin,
                limits.vUp,
                std::bind(&Parametrization::FunctionToDEdxIntegral, parametrization_, energy, std::placeholders::_1),
                4);
        }
    }

    return energy * sum;
}

// ------------------------------------------------------------------------- //
double EpairIntegral::FunctionToDEdxIntegralReverse(double energy, double v)
{
    return (1 - v) * parametrization_->DifferentialCrossSection(energy, v);
}
