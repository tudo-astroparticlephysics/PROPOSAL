
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
            
        sum += dedx_integral_.Integrate(
            limits.vMin,
            limits.vUp,
            std::bind(&Parametrization::FunctionToDEdxIntegral, parametrization_, energy, std::placeholders::_1),
            4);
    }

    return energy * sum;
}


