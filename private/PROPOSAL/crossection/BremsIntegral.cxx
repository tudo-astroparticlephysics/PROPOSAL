
#include <functional>

#include "PROPOSAL/crossection/BremsIntegral.h"
#include "PROPOSAL/crossection/parametrization/Bremsstrahlung.h"
#include "PROPOSAL/medium/Medium.h"

using namespace PROPOSAL;

BremsIntegral::BremsIntegral(const Bremsstrahlung& param, std::shared_ptr<EnergyCutSettings> cuts)
    : CrossSectionIntegral(param, cuts)
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
    double vUp;

    for (int i = 0; i < (parametrization_->GetMedium().GetNumComponents()); i++)
    {
        parametrization_->SetCurrentComponent(i);
        Parametrization::KinematicLimits limits = parametrization_->GetKinematicLimits(energy);
        vUp = cuts_.GetCut(energy);

        sum += dedx_integral_.Integrate(
            limits.vMin,
            vUp,
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
