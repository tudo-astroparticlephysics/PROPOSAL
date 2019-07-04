
#include <functional>
#include <cmath>

#include "PROPOSAL/crossection/CrossSectionInterpolant.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/math/Interpolant.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

CrossSectionInterpolant::CrossSectionInterpolant(const DynamicData::Type& type, const Parametrization& param)
    : CrossSection(type, param)
    , dedx_interpolant_(NULL)
    , de2dx_interpolant_(NULL)
    , dndx_interpolant_1d_(param.GetMedium().GetNumComponents(), NULL)
    , dndx_interpolant_2d_(param.GetMedium().GetNumComponents(), NULL)
{
}

bool CrossSectionInterpolant::compare(const CrossSection& cross_section) const
{
    const CrossSectionInterpolant* cross_section_interpolant =
        static_cast<const CrossSectionInterpolant*>(&cross_section);

    if (*dedx_interpolant_ != *cross_section_interpolant->dedx_interpolant_)
        return false;
    else if (*de2dx_interpolant_ != *cross_section_interpolant->de2dx_interpolant_)
        return false;
    else if (dndx_interpolant_1d_.size() != cross_section_interpolant->dndx_interpolant_1d_.size())
        return false;
    else if (dndx_interpolant_2d_.size() != cross_section_interpolant->dndx_interpolant_2d_.size())
        return false;

    for (unsigned int i = 0; i < dndx_interpolant_1d_.size(); ++i)
    {
        if (*dndx_interpolant_1d_[i] != *cross_section_interpolant->dndx_interpolant_1d_[i])
            return false;
    }
    for (unsigned int i = 0; i < dndx_interpolant_2d_.size(); ++i)
    {
        if (*dndx_interpolant_2d_[i] != *cross_section_interpolant->dndx_interpolant_2d_[i])
            return false;
    }

    return true;
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
        builder2d[i]
            .SetMax1(def.nodes_cross_section)
            .SetX1Min(parametrization_->GetParticleDef().mass)
            .SetX1Max(def.max_node_energy)
            .SetMax2(def.nodes_cross_section)
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
            .SetFunction2D(std::bind(
                &CrossSectionInterpolant::FunctionToBuildDNdxInterpolant2D,
                this,
                std::placeholders::_1,
                std::placeholders::_2,
                std::ref(integral),
                i));

        builder_container2d[i].first  = &builder2d[i];
        builder_container2d[i].second = &dndx_interpolant_2d_[i];

        builder1d[i]
            .SetMax(def.nodes_cross_section)
            .SetXMin(parametrization_->GetParticleDef().mass)
            .SetXMax(def.max_node_energy)
            .SetRomberg(def.order_of_interpolation)
            .SetRational(false)
            .SetRelative(false)
            .SetIsLog(true)
            .SetRombergY(def.order_of_interpolation)
            .SetRationalY(true)
            .SetRelativeY(false)
            .SetLogSubst(false)
            .SetFunction1D(std::bind(&CrossSectionInterpolant::FunctionToBuildDNdxInterpolant, this, std::placeholders::_1, i));

        builder_container1d[i].first  = &builder1d[i];
        builder_container1d[i].second = &dndx_interpolant_1d_[i];
    }

    builder_return.insert(builder_return.end(), builder_container2d.begin(), builder_container2d.end());
    builder_return.insert(builder_return.end(), builder_container1d.begin(), builder_container1d.end());
    // builder2d.insert(builder2d.end(), builder1d.begin(), builder1d.end());

    Helper::InitializeInterpolation("dNdx", builder_return, std::vector<Parametrization*>(1, parametrization_), def);
}

CrossSectionInterpolant::CrossSectionInterpolant(const CrossSectionInterpolant& cross_section)
    : CrossSection(cross_section)
{
    if (cross_section.dedx_interpolant_ != NULL)
    {
        dedx_interpolant_ = new Interpolant(*cross_section.dedx_interpolant_);
    }
    else{
        dedx_interpolant_ = NULL;
    }

    if (cross_section.de2dx_interpolant_ != NULL)
    {
        de2dx_interpolant_ = new Interpolant(*cross_section.de2dx_interpolant_);
    }
    else{
        de2dx_interpolant_ = NULL;
    }

    int num_components = cross_section.parametrization_->GetMedium().GetNumComponents();

    dndx_interpolant_1d_.reserve(num_components);
    for (auto interpolant: cross_section.dndx_interpolant_1d_)
    {
        if (interpolant != NULL)
        {
            dndx_interpolant_1d_.push_back(new Interpolant(*interpolant));
        }
    }

    dndx_interpolant_2d_.reserve(num_components);
    for (auto interpolant: cross_section.dndx_interpolant_2d_)
    {
        if (interpolant != NULL)
        {
            dndx_interpolant_2d_.push_back(new Interpolant(*interpolant));
        }
    }
}

CrossSectionInterpolant::~CrossSectionInterpolant()
{
    delete dedx_interpolant_;
    delete de2dx_interpolant_;

    for (auto interpolant: dndx_interpolant_1d_)
    {
        delete interpolant;
    }

    for (auto interpolant: dndx_interpolant_2d_)
    {
        delete interpolant;
    }
}

// ------------------------------------------------------------------------- //
// Pulblic methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double CrossSectionInterpolant::CalculatedE2dx(double energy)
{
    return parametrization_->GetMultiplier() * std::max(de2dx_interpolant_->Interpolate(energy), 0.0);
}

// ------------------------------------------------------------------------- //
double CrossSectionInterpolant::CalculatedNdx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    sum_of_rates_ = 0;

    const ComponentVec& components = parametrization_->GetMedium().GetComponents();
    for (size_t i = 0; i < components.size(); ++i)
    {
        prob_for_component_[i] = std::max(dndx_interpolant_1d_[i]->Interpolate(energy), 0.);
        sum_of_rates_ += prob_for_component_[i];
    }
    return parametrization_->GetMultiplier() * sum_of_rates_;
}

// ------------------------------------------------------------------------- //
double CrossSectionInterpolant::CalculatedNdx(double energy, double rnd)
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

    const ComponentVec& components = parametrization_->GetMedium().GetComponents();
    for (size_t i = 0; i < components.size(); ++i)
    {
        prob_for_component_[i] = std::max(dndx_interpolant_1d_[i]->Interpolate(energy), 0.);
        sum_of_rates_ += prob_for_component_[i];
    }

    return parametrization_->GetMultiplier() * sum_of_rates_;
}

// ------------------------------------------------------------------------- //
double CrossSectionInterpolant::CalculateStochasticLoss(double energy, double rnd1, double rnd2)
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
double CrossSectionInterpolant::CalculateStochasticLoss(double energy, double rnd1)
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
            Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);

            if (limits.vUp == limits.vMax)
            {
                return energy * limits.vUp;
            }

            return energy *
                   (limits.vUp * std::exp(dndx_interpolant_2d_.at(i)->FindLimit(energy, rnd_ * prob_for_component_[i]) *
                                     std::log(limits.vMax / limits.vUp)));
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
    return 0; // just to prevent warnings
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
double CrossSectionInterpolant::FunctionToBuildDNdxInterpolant2D(double energy,
                                                                 double v,
                                                                 Integral& integral,
                                                                 int component)
{
    parametrization_->SetCurrentComponent(component);
    Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);

    if (limits.vUp == limits.vMax)
    {
        return 0;
    }

    v = limits.vUp * std::exp(v * std::log(limits.vMax / limits.vUp));

    return integral.Integrate(
        limits.vUp, v, std::bind(&Parametrization::FunctionToDNdxIntegral, parametrization_, energy, std::placeholders::_1), 4);
}
