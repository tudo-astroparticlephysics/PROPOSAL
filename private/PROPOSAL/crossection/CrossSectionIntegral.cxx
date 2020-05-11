
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

CrossSectionIntegral::CrossSectionIntegral(unique_ptr<Parametrization>&& param,
    shared_ptr<const EnergyCutSettings> cut)
    : CrossSection(std::forward<unique_ptr<Parametrization>>(param), cut)
{
    for (auto comp : parametrization_->GetComponents()) {
        parametrization_->SetCurrentComponent(comp);
        dedx_integral_.emplace_back(
            bind(&CrossSectionIntegral::dedx_integral, this, _1));
        de2dx_integral_.emplace_back(
            bind(&CrossSectionIntegral::de2dx_integral, this, _1));
        dndx_integral_[comp.GetHash()]
            = bind(&CrossSectionIntegral::dndx_integral, this, _1, _2);
    }
}

double CrossSectionIntegral::dndx_integral(double energy, double rnd)
{
    auto v_cut = GetEnergyCut(energy);
    auto v_max = parametrization_->GetKinematicLimits(energy).vMax;

    auto rate = integral_.IntegrateWithRandomRatio(v_cut, v_max,
        bind(&Parametrization::FunctionToDNdxIntegral, parametrization_.get(),
            energy, _1),
        4, rnd);

    if (rnd != 1)
        return integral_.GetUpperLimit();

    return rate;
}

double CrossSectionIntegral::dedx_integral(double energy)
{
    auto v_cut = GetEnergyCut(energy);
    auto v_min = parametrization_->GetKinematicLimits(energy).vMin;

    return integral_.Integrate(v_min, v_cut,
        bind(&Parametrization::FunctionToDEdxIntegral, parametrization_.get(),
            energy, _1),
        4);
}

double CrossSectionIntegral::de2dx_integral(double energy)
{
    auto v_min = parametrization_->GetKinematicLimits(energy).vMin;
    auto v_cut = GetEnergyCut(energy);

    return integral_.Integrate(v_min, v_cut,
        bind(&Parametrization::FunctionToDE2dxIntegral, parametrization_.get(),
            energy, _1),
        2);
}

unordered_map<size_t, double> CrossSectionIntegral::CalculatedNdx(double energy)
{
    unordered_map<size_t, double> rates;
    for (auto& comp : parametrization_->GetComponents())
        rates[comp.GetHash()] = dndx_integral_[comp.GetHash()](energy, 1.);

    return rates;
}

double CrossSectionIntegral::CalculateStochasticLoss(
    double energy, double rate, size_t comp_hash)
{
    return dndx_integral_[comp_hash](energy, -rate);
}

double CrossSectionIntegral::CalculatedEdx(double energy)
{
    double sum = 0.;
    for (const auto& dedx : dedx_integral_)
        sum += dedx(energy);

    return energy * sum;
}

double CrossSectionIntegral::CalculatedE2dx(double energy)
{
    double sum = 0.;
    for (const auto& de2dx : de2dx_integral_)
        sum += de2dx(energy);

    return energy * energy * sum;
}
