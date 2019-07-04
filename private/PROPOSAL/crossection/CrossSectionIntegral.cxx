
#include <functional>

#include "PROPOSAL/crossection/CrossSectionIntegral.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

CrossSectionIntegral::CrossSectionIntegral(const DynamicData::Type& type, const Parametrization& param)
    : CrossSection(type, param)
    , dedx_integral_(IROMB, IMAXS, IPREC)
    , de2dx_integral_(IROMB, IMAXS, IPREC)
    , dndx_integral_(param.GetMedium().GetNumComponents(), Integral(IROMB, IMAXS, IPREC))
{
}

CrossSectionIntegral::CrossSectionIntegral(const CrossSectionIntegral& cross_section)
    : CrossSection(cross_section)
    , dedx_integral_(cross_section.dedx_integral_)
    , de2dx_integral_(cross_section.de2dx_integral_)
    , dndx_integral_(cross_section.dndx_integral_)
{
}

CrossSectionIntegral::~CrossSectionIntegral() {}

bool CrossSectionIntegral::compare(const CrossSection& cross_section) const
{
    const CrossSectionIntegral* cross_section_integral = static_cast<const CrossSectionIntegral*>(&cross_section);

    if (dedx_integral_ != cross_section_integral->dedx_integral_)
        return false;
    else if (de2dx_integral_ != cross_section_integral->de2dx_integral_)
        return false;
    else if (dndx_integral_ != cross_section_integral->dndx_integral_)
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
// Pulblic methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double CrossSectionIntegral::CalculatedE2dx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    double aux = 0;
    aux = CrossSectionIntegral::CalculatedE2dxWithoutMultiplier(energy);


    return parametrization_->GetMultiplier() * aux;
}

double CrossSectionIntegral::CalculatedE2dxWithoutMultiplier(double energy)
{
    double sum = 0;

    for (size_t i = 0; i < components_.size(); ++i)
    {
        parametrization_->SetCurrentComponent(i);
        Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);

        sum += de2dx_integral_.Integrate(
            limits.vMin,
            limits.vUp,
            std::bind(&Parametrization::FunctionToDE2dxIntegral, parametrization_, energy, std::placeholders::_1),
            2);
    }

    return energy * energy * sum;
}

// ------------------------------------------------------------------------- //
double CrossSectionIntegral::CalculatedNdx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    sum_of_rates_ = 0;

    for (size_t i = 0; i < components_.size(); ++i)
    {
        parametrization_->SetCurrentComponent(i);
        Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);

        prob_for_component_[i] = dndx_integral_[i].Integrate(
            limits.vUp,
            limits.vMax,
            std::bind(&Parametrization::FunctionToDNdxIntegral, parametrization_, energy, std::placeholders::_1),
            4);
        sum_of_rates_ += prob_for_component_[i];
    }
    return parametrization_->GetMultiplier() * sum_of_rates_;
}

// ------------------------------------------------------------------------- //
double CrossSectionIntegral::CalculatedNdx(double energy, double rnd)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    // The random number will be stored to be able
    // to check if dNdx is already calculated for this random number.
    // This avoids a second calculation in CalculateStochaticLoss
    rnd_ = rnd;

    sum_of_rates_ = 0;

    for (size_t i = 0; i < components_.size(); ++i)
    {
        parametrization_->SetCurrentComponent(i);
        Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);

        prob_for_component_[i] = dndx_integral_[i].IntegrateWithRandomRatio(
            limits.vUp,
            limits.vMax,
            std::bind(&Parametrization::FunctionToDNdxIntegral, parametrization_, energy, std::placeholders::_1),
            4,
            rnd);
        sum_of_rates_ += prob_for_component_[i];
    }

    return parametrization_->GetMultiplier() * sum_of_rates_;
}

// ------------------------------------------------------------------------- //
double CrossSectionIntegral::CalculateStochasticLoss(double energy, double rnd1, double rnd2)
{
    if (rnd1 != rnd_)
    {
        CalculatedNdx(energy, rnd1);
    }

    return CalculateStochasticLoss(energy, rnd2);
}

// ------------------------------------------------------------------------- //
// Private methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double CrossSectionIntegral::CalculateStochasticLoss(double energy, double rnd1)
{

    double rnd;
    double rsum;

    rnd  = rnd1 * sum_of_rates_;
    rsum = 0;

    for (size_t i = 0; i < components_.size(); ++i)
    {
        rsum += prob_for_component_[i];

        if (rsum > rnd)
        {
            parametrization_->SetCurrentComponent(i);
            return energy * dndx_integral_[i].GetUpperLimit();
        }
    }

    // sometime everything is fine, just the probability for interaction is zero
    bool prob_for_all_comp_is_zero = true;
    for (size_t i = 0; i < components_.size(); ++i)
    {
        parametrization_->SetCurrentComponent(i);
        Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);

        if (limits.vUp != limits.vMax)
            prob_for_all_comp_is_zero = false;
    }

    if (prob_for_all_comp_is_zero)
        return 0;

    log_fatal("sum was not initialized correctly");
    return 0;
}
