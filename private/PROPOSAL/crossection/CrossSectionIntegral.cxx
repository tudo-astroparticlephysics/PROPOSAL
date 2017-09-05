
#include <boost/bind.hpp>

#include "PROPOSAL/crossection/CrossSectionIntegral.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/Output.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

CrossSectionIntegral::CrossSectionIntegral(Parametrization& param)
    : CrossSection(param)
{
}

CrossSectionIntegral::CrossSectionIntegral(const CrossSectionIntegral& cross_section)
    : CrossSection(cross_section)
{
}

CrossSectionIntegral::~CrossSectionIntegral()
{
}

// ------------------------------------------------------------------------- //
// Pulblic methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double CrossSectionIntegral::CalculatedNdx(double energy)
{
    if(parametrization.GetMultiplier() <= 0)
    {
        return 0;
    }

    sum_of_rates_ = 0;

    const ComponentVec& components = parametrization.GetMedium().GetComponents();
    for(size_t i = 0; i < components.size(); ++i)
    {
        parametrization.SetCurrentComponent(components[i]);
        Parametrization::IntegralLimits limits = parametrization.GetIntegralLimits(energy);

        prob_for_component_[i] = dndx_integral_[i].Integrate(limits.vUp, limits.vMax, boost::bind(&Parametrization::FunctionToDNdxIntegral, &parametrization, energy,  _1),4);
        sum_of_rates_ += prob_for_component_[i];
    }
    return sum_of_rates_;
}

// ------------------------------------------------------------------------- //
double CrossSectionIntegral::CalculatedNdx(double energy, double rnd)
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
        parametrization.SetCurrentComponent(components[i]);
        Parametrization::IntegralLimits limits = parametrization.GetIntegralLimits(energy);

        prob_for_component_.at(i) = dndx_integral_[i].IntegrateWithRandomRatio(limits.vUp, limits.vMax, boost::bind(&Parametrization::FunctionToDNdxIntegral, &parametrization, energy,  _1), 4, rnd);
        sum_of_rates_ += prob_for_component_.at(i);
    }

    return sum_of_rates_;
}

// ------------------------------------------------------------------------- //
double CrossSectionIntegral::CalculateStochasticLoss(double energy,  double rnd1, double rnd2)
{
    if(rnd1 != rnd_ )
    {
        CalculatedNdx(energy, rnd1);
    }

    return CalculateStochasticLoss(energy, rnd2);
}

// ------------------------------------------------------------------------- //
double CrossSectionIntegral::FunctionToBuildDNdxInterpolant(double energy, int component)
{
    (void) energy;
    (void) component;
    // Empty, no interpolation needed
}

// ------------------------------------------------------------------------- //
double CrossSectionIntegral::FunctionToBuildDNdxInterpolant2D(double energy, double v, int component)
{
    (void) energy;
    (void) v;
    (void) component;
    // Empty, no interpolation needed
}

// ------------------------------------------------------------------------- //
// Private methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double CrossSectionIntegral::CalculateStochasticLoss(double energy, double rnd1)
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
            return energy * dndx_integral_[i].GetUpperLimit();
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
