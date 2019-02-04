
#include <functional>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/MupairIntegral.h"
#include "PROPOSAL/crossection/parametrization/MupairProduction.h"
#include "PROPOSAL/medium/Medium.h"

using namespace PROPOSAL;

MupairIntegral::MupairIntegral(const MupairProduction& param)
    : CrossSectionIntegral(DynamicData::MuPair, param)
{
}

MupairIntegral::MupairIntegral(const MupairIntegral& mupair)
    : CrossSectionIntegral(mupair)
{
}

MupairIntegral::~MupairIntegral() {}

// ----------------------------------------------------------------- //
// Public methods
// ----------------------------------------------------------------- //

double MupairIntegral::CalculatedEdx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

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
                                         std::bind(&MupairIntegral::FunctionToDEdxIntegralReverse, this, energy, std::placeholders::_1),
                                         2) +
                dedx_integral_.Integrate(
                    r2, 1 - r1, std::bind(&MupairIntegral::FunctionToDEdxIntegralReverse, this, energy, std::placeholders::_1), 4);

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
double MupairIntegral::FunctionToDEdxIntegralReverse(double energy, double v)
{
    return (1 - v) * parametrization_->DifferentialCrossSection(energy, v);
}
