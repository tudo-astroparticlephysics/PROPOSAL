
#include <boost/bind.hpp>
#include <cmath>

#include "PROPOSAL/crossection/CrossSectionInterpolant.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/math/Interpolant.h"

#include "PROPOSAL/Output.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

CrossSectionInterpolant::CrossSectionInterpolant(const DynamicData::Type& type, const Parametrization& param, InterpolationDef def)
    : CrossSection(type, param)
    , dedx_interpolant_(NULL)
    , de2dx_interpolant_(NULL)
    , dndx_interpolant_1d_(param.GetMedium().GetNumComponents(), NULL)
    , dndx_interpolant_2d_(param.GetMedium().GetNumComponents(), NULL)
{
    // InitdNdxInerpolation(def);
}

// ------------------------------------------------------------------------- //
void CrossSectionInterpolant::InitdNdxInerpolation(const InterpolationDef& def)
{
    // --------------------------------------------------------------------- //
    // Builder for dNdx
    // --------------------------------------------------------------------- //

    std::vector<Interpolant1DBuilder> builder1d(components_.size());
    std::vector<Interpolant2DBuilder> builder2d(components_.size());

    Helper::InterpolantBuilderContainer builder_container1d(components_.size());
    Helper::InterpolantBuilderContainer builder_container2d(components_.size());
    Helper::InterpolantBuilderContainer builder_return;

    Integral integral(IROMB, IMAXS, IPREC);

    for (unsigned int i = 0; i < components_.size(); ++i)
    {
        // !!! IMPORTANT !!!
        // Order of builder matter because the functions needed for 1d interpolation
        // needs the already intitialized 2d interpolants.
        builder2d[i].SetMax1(NUM1)
            .SetX1Min(parametrization_->GetParticleDef().low)
            .SetX1Max(BIGENERGY)
            .SetMax2(NUM1)
            .SetX2Min(0.0)
            .SetX2Max(1.0)
            .SetRomberg1(def.order_of_interpolation)
            .SetRational1(false)
            .SetRelative1(false)
            .SetIsLog1(true)
            .SetRomberg2(def.order_of_interpolation)
            .SetRational2(false)
            .SetRelative2(false)
            .SetIsLog2(false)
            .SetRombergY(def.order_of_interpolation)
            .SetRationalY(true)
            .SetRelativeY(false)
            .SetLogSubst(false)
            .SetFunction2D(boost::bind(&CrossSectionInterpolant::FunctionToBuildDNdxInterpolant2D, this, _1, _2, boost::ref(integral), i));

        builder_container2d[i].first = &builder2d[i];
        builder_container2d[i].second = &dndx_interpolant_2d_[i];

        builder1d[i].SetMax(NUM1)
            .SetXMin(parametrization_->GetParticleDef().low)
            .SetXMax(BIGENERGY)
            .SetRomberg(def.order_of_interpolation)
            .SetRational(false)
            .SetRelative(false)
            .SetIsLog(true)
            .SetRombergY(def.order_of_interpolation)
            .SetRationalY(true)
            .SetRelativeY(false)
            .SetLogSubst(false)
            .SetFunction1D(boost::bind(&CrossSectionInterpolant::FunctionToBuildDNdxInterpolant, this, _1, i));

        builder_container1d[i].first = &builder1d[i];
        builder_container1d[i].second = &dndx_interpolant_1d_[i];
    }

    builder_return.insert(builder_return.end(), builder_container2d.begin(), builder_container2d.end());
    builder_return.insert(builder_return.end(), builder_container1d.begin(), builder_container1d.end());
    // builder2d.insert(builder2d.end(), builder1d.begin(), builder1d.end());

    Helper::InitializeInterpolation("dNdx",
                                    builder_return,
                                    std::vector<Parametrization*>(1, parametrization_), def);
}

CrossSectionInterpolant::CrossSectionInterpolant(const CrossSectionInterpolant& cross_section)
    : CrossSection(cross_section)
{
    if (cross_section.dedx_interpolant_ != NULL)
    {
        dedx_interpolant_ = new Interpolant(*cross_section.dedx_interpolant_);
    }

    if (cross_section.de2dx_interpolant_ != NULL)
    {
        de2dx_interpolant_ = new Interpolant(*cross_section.de2dx_interpolant_);
    }

    int num_components = cross_section.parametrization_->GetMedium().GetNumComponents();

    dndx_interpolant_1d_.reserve(num_components);
    for(InterpolantVec::const_iterator iter = cross_section.dndx_interpolant_1d_.begin(); iter != cross_section.dndx_interpolant_1d_.end(); ++iter)
    {
        if (*iter != NULL)
        {
            dndx_interpolant_1d_.push_back(new Interpolant(**iter));
        }
    }

    dndx_interpolant_2d_.reserve(num_components);
    for(InterpolantVec::const_iterator iter = cross_section.dndx_interpolant_2d_.begin(); iter != cross_section.dndx_interpolant_2d_.end(); ++iter)
    {
        if (*iter != NULL)
        {
            dndx_interpolant_2d_.push_back(new Interpolant(**iter));
        }
    }
}

CrossSectionInterpolant::~CrossSectionInterpolant()
{
    delete dedx_interpolant_;
    delete de2dx_interpolant_;

    for (InterpolantVec::const_iterator iter = dndx_interpolant_1d_.begin(); iter != dndx_interpolant_1d_.end(); ++iter) {
        delete *iter;
    }

    for (InterpolantVec::const_iterator iter = dndx_interpolant_2d_.begin(); iter != dndx_interpolant_2d_.end(); ++iter) {
        delete *iter;
    }
}

// ------------------------------------------------------------------------- //
// Pulblic methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double CrossSectionInterpolant::CalculatedE2dx(double energy)
{
    return std::max(de2dx_interpolant_->Interpolate(energy), 0.0);
}

// ------------------------------------------------------------------------- //
double CrossSectionInterpolant::CalculatedNdx(double energy)
{
    if(parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    sum_of_rates_ = 0;

    const ComponentVec& components = parametrization_->GetMedium().GetComponents();
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
    if(parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    // The random number will be stored to be able
    // to check if dNdx is already calculated for this random number.
    // This avoids a second calculation in CalculateStochaticLoss
    rnd_ = rnd;

    sum_of_rates_ = 0;

    const ComponentVec& components = parametrization_->GetMedium().GetComponents();
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

    for(size_t i = 0; i < components_.size(); ++i)
    {
        rsum    += prob_for_component_[i];

        if(rsum > rnd)
        {
            parametrization_->SetCurrentComponent(i);
            Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);

            if(limits.vUp == limits.vMax)
            {
                return energy * limits.vUp;
            }

            return energy * (limits.vUp * exp(dndx_interpolant_2d_.at(i)->FindLimit(energy, rnd_ * prob_for_component_[i])*log(limits.vMax / limits.vUp)));
        }
    }

    //sometime everything is fine, just the probability for interaction is zero
    bool prob_for_all_comp_is_zero=true;
    for(size_t i = 0; i < components_.size(); ++i)
    {
        parametrization_->SetCurrentComponent(i);
        Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);

        if (limits.vUp != limits.vMax)
            prob_for_all_comp_is_zero = false;
    }

    if (prob_for_all_comp_is_zero)
        return 0;

    log_fatal("sum was not initialized correctly");
}

// ------------------------------------------------------------------------- //
// Function needed for interpolation intitialization
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double CrossSectionInterpolant::FunctionToBuildDNdxInterpolant(double energy, int component)
{
    return dndx_interpolant_2d_[component]->Interpolate(energy, 1.);
}

//----------------------------------------------------------------------------//
double CrossSectionInterpolant::FunctionToBuildDNdxInterpolant2D(double energy, double v, Integral& integral, int component)
{
    parametrization_->SetCurrentComponent(component);
    Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);

    if (limits.vUp == limits.vMax)
    {
        return 0;
    }

    v = limits.vUp * exp(v * log(limits.vMax / limits.vUp));

    return integral.Integrate(
        limits.vUp, v, boost::bind(&Parametrization::FunctionToDNdxIntegral, parametrization_, energy, _1), 4);
}
