
#include <functional>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/PhotoIntegral.h"
#include "PROPOSAL/crossection/parametrization/Photonuclear.h"
#include "PROPOSAL/medium/Medium.h"

using namespace PROPOSAL;

PhotoIntegral::PhotoIntegral(const Photonuclear& param)
    : CrossSectionIntegral(DynamicData::NuclInt, param)
{
}

PhotoIntegral::PhotoIntegral(const PhotoIntegral& brems)
    : CrossSectionIntegral(brems)
{
}

PhotoIntegral::~PhotoIntegral() {}

// ----------------------------------------------------------------- //
// Public methods
// ----------------------------------------------------------------- //

double PhotoIntegral::CalculatedEdx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    return parametrization_->GetMultiplier() * PhotoIntegral::CalculatedEdxWithoutMultiplier(energy);
}

double PhotoIntegral::CalculatedEdxWithoutMultiplier(double energy)
{
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
