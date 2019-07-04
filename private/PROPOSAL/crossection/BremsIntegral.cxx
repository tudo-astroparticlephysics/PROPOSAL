
#include <functional>

#include "PROPOSAL/crossection/BremsIntegral.h"
#include "PROPOSAL/crossection/parametrization/Bremsstrahlung.h"
#include "PROPOSAL/medium/Medium.h"

using namespace PROPOSAL;

BremsIntegral::BremsIntegral(const Bremsstrahlung& param)
    : CrossSectionIntegral(DynamicData::Brems, param)
{
}

BremsIntegral::BremsIntegral(const BremsIntegral& brems)
    : CrossSectionIntegral(brems)
{
}

BremsIntegral::~BremsIntegral() {}

// ----------------------------------------------------------------- //
// Public methods
// ----------------------------------------------------------------- //
double BremsIntegral::CalculatedEdxWithoutMultiplier(double energy){
    double sum = 0;

    for (int i = 0; i < (parametrization_->GetMedium().GetNumComponents()); i++)
    {
        parametrization_->SetCurrentComponent(i);
        Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);

        sum += dedx_integral_.Integrate(
            limits.vMin,
            limits.vUp,
            std::bind(&Parametrization::FunctionToDEdxIntegral, parametrization_, energy, std::placeholders::_1),
            2);
    }

    return energy * sum;
}

double BremsIntegral::CalculatedEdx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    return parametrization_->GetMultiplier() * BremsIntegral::CalculatedEdxWithoutMultiplier(energy);
}
