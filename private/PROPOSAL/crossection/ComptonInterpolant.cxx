
#include <functional>
#include <cmath>

#include "PROPOSAL/crossection/ComptonIntegral.h"
#include "PROPOSAL/crossection/ComptonInterpolant.h"
#include "PROPOSAL/crossection/parametrization/Compton.h"

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/math/RandomGenerator.h"

using namespace PROPOSAL;

ComptonInterpolant::ComptonInterpolant(const Compton& param, InterpolationDef def)
        : CrossSectionInterpolant(InteractionType::Compton, param)
{
    // Use own CrossSecition dNdx interpolation
    ComptonInterpolant::InitdNdxInterpolation(def);

    // --------------------------------------------------------------------- //
    // Builder for DEdx
    // --------------------------------------------------------------------- //

    Interpolant1DBuilder builder1d;
    Helper::InterpolantBuilderContainer builder_container;

    // Needed for CalculatedEdx integration
    ComptonIntegral compton(param);

    builder1d.SetMax(def.nodes_cross_section)
            .SetXMin(param.GetParticleDef().low)
            .SetXMax(def.max_node_energy)
            .SetRomberg(def.order_of_interpolation)
            .SetRational(true)
            .SetRelative(false)
            .SetIsLog(true)
            .SetRombergY(def.order_of_interpolation)
            .SetRationalY(false)
            .SetRelativeY(false)
            .SetLogSubst(true)
            .SetFunction1D(std::bind(&CrossSectionIntegral::CalculatedEdxWithoutMultiplier, &compton, std::placeholders::_1));

    builder_container.push_back(std::make_pair(&builder1d, &dedx_interpolant_));

    // --------------------------------------------------------------------- //
    // Builder for DE2dx
    // --------------------------------------------------------------------- //

    Interpolant1DBuilder builder_de2dx;
    Helper::InterpolantBuilderContainer builder_container_de2dx;

    builder_de2dx.SetMax(def.nodes_continous_randomization)
            .SetXMin(param.GetParticleDef().low)
            .SetXMax(def.max_node_energy)
            .SetRomberg(def.order_of_interpolation)
            .SetRational(false)
            .SetRelative(false)
            .SetIsLog(true)
            .SetRombergY(def.order_of_interpolation)
            .SetRationalY(false)
            .SetRelativeY(false)
            .SetLogSubst(false)
            .SetFunction1D(std::bind(&CrossSectionIntegral::CalculatedE2dxWithoutMultiplier, &compton, std::placeholders::_1));

    builder_container_de2dx.push_back(std::make_pair(&builder_de2dx, &de2dx_interpolant_));

    Helper::InitializeInterpolation("dEdx", builder_container, std::vector<Parametrization*>(1, parametrization_), def);
    Helper::InitializeInterpolation(
            "dE2dx", builder_container_de2dx, std::vector<Parametrization*>(1, parametrization_), def);
}

ComptonInterpolant::ComptonInterpolant(const ComptonInterpolant& compton)
        : CrossSectionInterpolant(compton)
{
}

ComptonInterpolant::~ComptonInterpolant() {}

// ----------------------------------------------------------------- //
// Public methods
// ----------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double ComptonInterpolant::CalculatedEdx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    return parametrization_->GetMultiplier() * std::max(dedx_interpolant_->Interpolate(energy), 0.0);
}

//----------------------------------------------------------------------------//
double ComptonInterpolant::FunctionToBuildDNdxInterpolant2D(double energy,
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

    v = limits.vUp + (limits.vMax - limits.vUp) * v;

    // Integrate with the substitution t = ln(1-v) to avoid numerical problems
    auto integrand_substitution = [&](double energy, double t){
        return std::exp(t) * parametrization_->FunctionToDNdxIntegral(energy, 1 - std::exp(t));
    };

    double t_min = std::log(1. - v);
    double t_max = std::log(1. - limits.vUp);


    return integral.Integrate(
            t_min,
            t_max,
            std::bind(integrand_substitution, energy, std::placeholders::_1),
            2);
}

double ComptonInterpolant::CalculateCumulativeCrossSection(double energy, int component, double v)
{
    parametrization_->SetCurrentComponent(component);
    Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);

    v = (v - limits.vUp) / (limits.vMax - limits.vUp);

    return dndx_interpolant_2d_.at(component)->Interpolate(energy, v);
}

std::pair<double, double> ComptonInterpolant::StochasticDeflection(double energy, double energy_loss) {
    double theta_deflect = RandomGenerator::Get().RandomDouble() * 2 * PI; // random azimuth
    double cosphi = 1. - (ME * (1. / (energy - energy_loss) - 1. / energy));

    return std::make_pair(cosphi, theta_deflect);
}

// ------------------------------------------------------------------------- //
// Private methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double ComptonInterpolant::CalculateStochasticLoss(double energy, double rnd1)
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

            // Linear interpolation in v
            return energy * ( limits.vUp + (limits.vMax - limits.vMin) * dndx_interpolant_2d_.at(i)->FindLimit(energy, rnd_ * prob_for_component_[i]) );
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
void ComptonInterpolant::InitdNdxInterpolation(const InterpolationDef& def)
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
                .SetX1Min(parametrization_->GetParticleDef().low)
                .SetX1Max(def.max_node_energy)
                .SetMax2(def.nodes_cross_section)
                .SetX2Min(1. / (2. * (1. - def.nodes_cross_section)))
                .SetX2Max((1. - 2. * def.nodes_cross_section) / (2. * (1. - def.nodes_cross_section)))
                .SetRomberg1(def.order_of_interpolation)
                .SetRational1(false)
                .SetRelative1(false)
                .SetIsLog1(true)
                .SetRomberg2(def.order_of_interpolation)
                .SetRational2(true)
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
                .SetXMin(parametrization_->GetParticleDef().low)
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
