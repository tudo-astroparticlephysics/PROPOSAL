
#include <boost/bind.hpp>
#include <cmath>

#include "PROPOSAL/crossection/CrossSectionInterpolant.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/Output.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

CrossSectionInterpolant::CrossSectionInterpolant(Parametrization& param)
    : CrossSection(param)
    , dedx_interpolant_(NULL)
    , dndx_interpolant_1d_(param.GetMedium().GetNumComponents(), NULL)
    , dndx_interpolant_2d_(param.GetMedium().GetNumComponents(), NULL)
{
}

CrossSectionInterpolant::CrossSectionInterpolant(const CrossSectionInterpolant& cross_section)
    : CrossSection(cross_section)
{
    if (cross_section.dedx_interpolant_ != NULL)
        dedx_interpolant_ = new Interpolant(*cross_section.dedx_interpolant_);

    dndx_interpolant_1d_.resize(cross_section.parametrization.GetMedium().GetNumComponents());
    for(InterpolantVec::iterator iter = dndx_interpolant_1d_.begin(); iter != dndx_interpolant_1d_.end(); ++iter)
    {
        if (*iter != NULL)
            dndx_interpolant_1d_.push_back(new Interpolant(**iter));
    }

    dndx_interpolant_2d_.resize(cross_section.parametrization.GetMedium().GetNumComponents());
    for(InterpolantVec::iterator iter = dndx_interpolant_2d_.begin(); iter != dndx_interpolant_2d_.end(); ++iter)
    {
        if (*iter != NULL)
            dndx_interpolant_2d_.push_back(new Interpolant(**iter));
    }
}

CrossSectionInterpolant::~CrossSectionInterpolant()
{
}

// ------------------------------------------------------------------------- //
// Pulblic methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double CrossSectionInterpolant::CalculatedNdx(double energy)
{
    if(parametrization.GetMultiplier() <= 0)
    {
        return 0;
    }

    sum_of_rates_ = 0;

    const ComponentVec& components = parametrization.GetMedium().GetComponents();
    for(size_t i = 0; i < components.size(); ++i)
    {
        prob_for_component_[i] = std::max(dndx_interpolant_1d_.at(i)->Interpolate(energy), 0.);
        sum_of_rates_ += prob_for_component_[i];
    }
    return sum_of_rates_;
}

// ------------------------------------------------------------------------- //
double CrossSectionInterpolant::CalculatedNdx(double energy, double rnd)
{
    if(parametrization.GetMultiplier() <= 0)
    {
        return 0;
    }

    // The random number will be stored to be able
    // to check if dNdx is already calculated for this random number.
    // This avoids a second calculation in CalculateStochaticLoss
    rnd_ = rnd;

    sum_of_rates_ = 0;

    const ComponentVec& components = parametrization.GetMedium().GetComponents();
    for(size_t i = 0; i < components.size(); ++i)
    {
        prob_for_component_[i] = std::max(dndx_interpolant_1d_.at(i)->Interpolate(energy), 0.);
        sum_of_rates_ += prob_for_component_.at(i);
    }

    return sum_of_rates_;
}

// ------------------------------------------------------------------------- //
double CrossSectionInterpolant::CalculateStochasticLoss(double energy,  double rnd1, double rnd2)
{
    if(rnd1 != rnd_ )
    {
        CalculatedNdx(energy, rnd1);
    }

    return CalculateStochasticLoss(energy, rnd2);
}


// ------------------------------------------------------------------------- //
// Private methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double CrossSectionInterpolant::CalculateStochasticLoss(double energy, double rnd1)
{

    double rnd;
    double rsum;

    rnd    =   rnd1*sum_of_rates_;
    rsum    =   0;

    const ComponentVec& components = parametrization.GetMedium().GetComponents();
    for(size_t i = 0; i < components.size(); ++i)
    {
        rsum    += prob_for_component_[i];

        if(rsum > rnd)
        {
            parametrization.SetCurrentComponent(components[i]);
            Parametrization::IntegralLimits limits = parametrization.GetIntegralLimits(energy);

            if(limits.vUp == limits.vMax)
            {
                return energy * limits.vUp;
            }

            return energy * (limits.vUp * exp(dndx_interpolant_2d_.at(i)->FindLimit(energy, rnd_ * prob_for_component_[i])*log(limits.vMax / limits.vUp)));
        }
    }

    //sometime everything is fine, just the probability for interaction is zero
    bool prob_for_all_comp_is_zero=true;
    for(size_t i = 0; i < components.size(); ++i)
    {
        parametrization.SetCurrentComponent(components[i]);
        Parametrization::IntegralLimits limits = parametrization.GetIntegralLimits(energy);

        if (limits.vUp != limits.vMax)
            prob_for_all_comp_is_zero = false;
    }

    if (prob_for_all_comp_is_zero)
        return 0;

    log_fatal("sum was not initialized correctly");
    return 0;
}
