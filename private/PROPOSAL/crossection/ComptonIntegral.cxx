
#include <functional>

#include "PROPOSAL/crossection/ComptonIntegral.h"
#include "PROPOSAL/crossection/parametrization/Compton.h"
#include "PROPOSAL/medium/Medium.h"

using namespace PROPOSAL;

ComptonIntegral::ComptonIntegral(const Compton& param)
        : CrossSectionIntegral(DynamicData::Compton, param)
{
}

ComptonIntegral::ComptonIntegral(const ComptonIntegral& compton)
        : CrossSectionIntegral(compton)
{
}

ComptonIntegral::~ComptonIntegral() {}

// ----------------------------------------------------------------- //
// Public methods
// ----------------------------------------------------------------- //
double ComptonIntegral::CalculatedEdxWithoutMultiplier(double energy){
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

double ComptonIntegral::CalculatedEdx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    return parametrization_->GetMultiplier() * ComptonIntegral::CalculatedEdxWithoutMultiplier(energy);
}
