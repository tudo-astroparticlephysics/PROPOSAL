
#include <functional>

#include "PROPOSAL/crossection/CrossSectionIntegral.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"

using std::bind;
using std::vector;
using std::placeholders::_1;

using namespace PROPOSAL;

double CrossSectionIntegral::dndx_integral(double energy, double rnd)
{
    auto v_cut = GetEnergyCut(energy);
    auto v_max = parametrization_->GetKinematicLimits(energy).vMax;
    /* auto dndx_func = bind(&Parametrization::FunctionToDNdxIntegral, */
    /*     parametrization_, energy, _1); */
    auto dndx_func = [&](double v){ return parametrization_->FunctionToDNdxIntegral(energy, v);};

    integral_.IntegrateWithRandomRatio(v_cut, v_max, dndx_func, 4, rnd);

    return integral_.GetUpperLimit();
}

double CrossSectionIntegral::dedx_integral(double energy)
{
    auto v_cut = GetEnergyCut(energy);
    auto v_min = parametrization_->GetKinematicLimits(energy).vMin;
    auto dedx_func = [&](double v){ return parametrization_->FunctionToDEdxIntegral(energy, v);};

    return integral_.Integrate(v_min, v_cut, dedx_func, 2);
}

double CrossSectionIntegral::de2dx_integral(double energy)
{
    auto v_cut = GetEnergyCut(energy);
    auto v_min = parametrization_->GetKinematicLimits(energy).vMin;
    auto de2dx_func = [&](double v){ return parametrization_->FunctionToDE2dxIntegral(energy, v);};

    return integral_.Integrate(v_min, v_cut, de2dx_func, 2);
}

double CrossSectionIntegral::CalculatedNdx(double energy, double rnd)
{
    vector<double> rates;
    for (auto& dndx : dndx_integral_) rates.push_back(dndx(energy, 1.));

    auto total_rate = accumulate(rates.begin(), rates.end(), 0);

    size_t nth_component = 0;
    for (const auto& rate : rates) {
            total_rate -= rate / rnd;
            if (total_rate < 0)
                return dndx_integral_.at(nth_component)(energy, total_rate);
            nth_component += 1;
    }

    return total_rate;
}

double CrossSectionIntegral::CalculatedEdx(double energy)
{
    auto integrate_and_sum = [energy](double sum, function<double(double)> func) {
        return sum + func(energy);
    };

    auto sum = accumulate(
        dedx_integral_.begin(), dedx_integral_.end(), 0, integrate_and_sum);

    return energy * sum;
}

double CrossSectionIntegral::CalculatedE2dx(double energy)
{
    auto integrate_and_sum = [energy](double sum, function<double(double)> func) {
        return sum + func(energy);
    };

    auto sum = accumulate(
        de2dx_integral_.begin(), de2dx_integral_.end(), 0, integrate_and_sum);

    return energy * energy * sum;
}
