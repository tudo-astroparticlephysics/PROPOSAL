
#include <cmath>
#include <functional>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/ComptonIntegral.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"

using namespace PROPOSAL;
using std::exp;
using std::log;
using std::vector;

// Integrate with the substitution t = ln(1-v) to avoid numerical problems
double ComptonIntegral::log_substitution(double v) const { return log(1 - v); }

double ComptonIntegral::dndx_integral(double energy, double rnd)
{
    auto v_cut = GetEnergyCut(energy);
    auto v_max = parametrization_->GetKinematicLimits(energy).vMax;

    auto t_min = log_substitution(v_cut);
    auto t_max = log_substitution(v_max);

    /* auto dndx_func = bind(&Compton::FunctionToDNdxIntegral, parametrization_,
     * energy, _1); */
    auto func_transformed = [&](double t) {
        return exp(t)
            * parametrization_->FunctionToDNdxIntegral(energy, 1 - exp(t));
    };

    integral_.IntegrateWithRandomRatio(t_min, t_max, func_transformed, 4, rnd);

    return integral_.GetUpperLimit();
}

double ComptonIntegral::dedx_integral(double energy)
{
    auto v_cut = GetEnergyCut(energy);
    auto v_max = parametrization_->GetKinematicLimits(energy).vMax;

    auto t_min = log_substitution(v_cut);
    auto t_max = log_substitution(v_max);

    /* auto dedx_func = bind( */
    /*     &Parametrization::FunctionToDEdxIntegral, parametrization_, energy, _1); */
    auto func_transformed
        = [&](double t) { return exp(t) * (parametrization_->FunctionToDEdxIntegral)(energy, 1 - exp(t)); };

    return integral_.Integrate(t_min, t_max, func_transformed, 2);
}

double ComptonIntegral::de2dx_integral(double energy)
{
    auto v_cut = GetEnergyCut(energy);
    auto v_max = parametrization_->GetKinematicLimits(energy).vMax;

    auto t_min = log_substitution(v_cut);
    auto t_max = log_substitution(v_max);

    /* auto de2dx_func = bind(&Parametrization::FunctionToDE2dxIntegral, */
    /*     parametrization_, energy, _1); */
    auto func_transformed = [&](double t) { return exp(t) * (parametrization_->FunctionToDE2dxIntegral)(energy, 1 - exp(t)); };

    return integral_.Integrate(t_min, t_max, func_transformed, 2);
}
