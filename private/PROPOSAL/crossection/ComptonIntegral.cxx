
#include <functional>
#include <cmath>

#include "PROPOSAL/crossection/ComptonIntegral.h"
#include "PROPOSAL/crossection/parametrization/Compton.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/RandomGenerator.h"

using namespace PROPOSAL;

ComptonIntegral::ComptonIntegral(const Compton& param)
        : CrossSectionIntegral(InteractionType::Compton, param)
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

    for (int i = 0; i < (parametrization_->GetMedium()->GetNumComponents()); i++)
    {
        parametrization_->SetCurrentComponent(i);
        Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);

        // Integrate with the substitution t = ln(1-v) to avoid numerical problems
        auto integrand_substitution = [&](double energy, double t){
            return std::exp(t) * parametrization_->FunctionToDEdxIntegral(energy, 1 - std::exp(t));
        };

        double t_min = std::log(1. - limits.vUp);
        double t_max = std::log(1. - limits.vMin);

        sum += dedx_integral_.Integrate(
                t_min,
                t_max,
                std::bind(integrand_substitution, energy, std::placeholders::_1),
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

double ComptonIntegral::CalculatedE2dxWithoutMultiplier(double energy)
{
    double sum = 0;

    for (size_t i = 0; i < components_.size(); ++i)
    {
        parametrization_->SetCurrentComponent(i);
        Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);

        // Integrate with the substitution t = ln(1-v) to avoid numerical problems
        auto integrand_substitution = [&](double energy, double t){
            return std::exp(t) * parametrization_->FunctionToDE2dxIntegral(energy, 1 - std::exp(t));
        };

        double t_min = std::log(1. - limits.vUp);
        double t_max = std::log(1. - limits.vMin);

        sum += de2dx_integral_.Integrate(
                t_min,
                t_max,
                std::bind(integrand_substitution, energy, std::placeholders::_1),
                2);
    }

    return energy * energy * sum;
}

// ------------------------------------------------------------------------- //
double ComptonIntegral::CalculatedNdx(double energy)
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

        // Integrate with the substitution t = ln(1-v) to avoid numerical problems
        auto integrand_substitution = [&](double energy, double t){
            return std::exp(t) * parametrization_->FunctionToDNdxIntegral(energy, 1 - std::exp(t));
        };

        double t_min = std::log(1. - limits.vMax);
        double t_max = std::log(1. - limits.vUp);

        prob_for_component_[i] = dndx_integral_[i].Integrate(
                t_min,
                t_max,
                std::bind(integrand_substitution, energy, std::placeholders::_1),
                2);

        sum_of_rates_ += prob_for_component_[i];
    }
    return parametrization_->GetMultiplier() * sum_of_rates_;
}

// ------------------------------------------------------------------------- //
double ComptonIntegral::CalculatedNdx(double energy, double rnd)
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

        // Integrate with the substitution t = ln(1-v) to avoid numerical problems
        // Has to be considered when evaluating the UpperLimit of the integral!

        auto integrand_substitution = [&](double energy, double t){
            return std::exp(t) * parametrization_->FunctionToDNdxIntegral(energy, 1 - std::exp(t));
        };

        double t_min = std::log(1. - limits.vMax);
        double t_max = std::log(1. - limits.vUp);

        // Switch limits to be able to evaluate the UpperLimit with the given substitution
        prob_for_component_[i] = -dndx_integral_[i].IntegrateWithRandomRatio(
                t_max,
                t_min,
                std::bind(integrand_substitution, energy, std::placeholders::_1),
                3,
                rnd);

        sum_of_rates_ += prob_for_component_[i];
    }

    return parametrization_->GetMultiplier() * sum_of_rates_;
}

double ComptonIntegral::CalculateCumulativeCrossSection(double energy, int i, double v)
{
    parametrization_->SetCurrentComponent(i);
    Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);

    auto integrand_substitution = [&](double energy, double t){
        return std::exp(t) * parametrization_->FunctionToDNdxIntegral(energy, 1 - std::exp(t));
    };

    double t_min = std::log(1. - v);
    double t_max = std::log(1. - limits.vUp);

    return dndx_integral_.at(i).Integrate(
            t_min,
            t_max,
            std::bind(integrand_substitution, energy, std::placeholders::_1),
            2);
}

std::pair<double, double> ComptonIntegral::StochasticDeflection(double energy, double energy_loss) {
    double theta_deflect = RandomGenerator::Get().RandomDouble() * 2 * PI; // random azimuth
    double cosphi = 1. - (ME * (1. / (energy - energy_loss) - 1. / energy));

    return std::make_pair(cosphi, theta_deflect);
}
// ------------------------------------------------------------------------- //
// Private methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double ComptonIntegral::CalculateStochasticLoss(double energy, double rnd1)
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
            // Resubstitution, v = 1 - exp(t)
            double upper_limit_in_v = 1. - std::exp(dndx_integral_[i].GetUpperLimit());
            return energy * upper_limit_in_v;
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
