
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
        dndx_integral_.emplace_back(
            bind(&CrossSectionIntegral::dndx_integral, this, _1, _2));
        dedx_integral_.emplace_back(
            bind(&CrossSectionIntegral::dedx_integral, this, _1));
        de2dx_integral_.emplace_back(
            bind(&CrossSectionIntegral::de2dx_integral, this, _1));
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
        2);
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

vector<double> CrossSectionIntegral::CalculatedNdx(double energy)
{
    vector<double> rates;
    for (auto& dndx : dndx_integral_)
        rates.push_back(dndx(energy, 1.));

    return rates;
}

vector<double> CrossSectionIntegral::CalculateEnergyLoss(
    double energy, const vector<double>& component_rates)
{
    vector<double> relativ_loss;
    auto component_rate = component_rates.cbegin();
    for (auto& dndx : dndx_integral_) {
        relativ_loss.push_back(dndx(energy, *component_rate));
        ++component_rate;
    }

    return relativ_loss;
}

double CrossSectionIntegral::CalculatedEdx(double energy)
{
    double sum = 0.;
    for (auto dedx: dedx_integral_)
        sum += dedx(energy);

    return energy * sum;
}

double CrossSectionIntegral::CalculatedE2dx(double energy)
{
    double sum = 0.;
    for (auto de2dx : de2dx_integral_)
        sum += de2dx(energy);

    return energy * energy * sum;
}
