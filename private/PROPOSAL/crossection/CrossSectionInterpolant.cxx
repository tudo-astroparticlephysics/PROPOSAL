
#include <cmath>
#include <functional>
#include <memory>

#include "PROPOSAL/crossection/CrossSectionInterpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/medium/Medium.h"

using namespace PROPOSAL;
using std::bind;
using std::placeholders::_1;
using std::placeholders::_2;

CrossSectionInterpolant::CrossSectionInterpolant(
    unique_ptr<Parametrization>&& param,
    shared_ptr<const EnergyCutSettings> cut, const InterpolationDef& def)
    : CrossSectionIntegral(forward<unique_ptr<Parametrization>>(param), cut)
    , hash_interpol_def(def.GetHash())
    , dedx_interpolant_(init_dedx_interpolation(def))
    , de2dx_interpolant_(init_de2dx_interpolation(def))
    , dndx_interpolants_(init_dndx_interpolation(def))
{
}

double CrossSectionInterpolant::transform_relativ_loss(
    double v, double energy) const
{
    auto v_cut = GetEnergyCut(energy);
    auto v_max = parametrization_->GetKinematicLimits(energy).vMax;

    return v_cut * std::exp(v * std::log(v_max / v_cut));
}

double CrossSectionInterpolant::dndx_integral(double energy, double v)
{
    auto v_cut = GetEnergyCut(energy);
    v = transform_relativ_loss(v, energy);

    return integral_.Integrate(v_cut, v,
        bind(&Parametrization::FunctionToDNdxIntegral, parametrization_.get(),
            energy, _1),
        4);
}

unique_ptr<Interpolant> CrossSectionInterpolant::init_dedx_interpolation(
    const InterpolationDef& def)
{
    Interpolant1DBuilder::Definition interpol_def;
    interpol_def.function1d = [this](double energy) {
        return CrossSectionIntegral::CalculatedEdx(energy);
    };
    interpol_def.max = def.nodes_cross_section;
    interpol_def.xmin = parametrization_->GetLowerEnergyLim();
    interpol_def.xmax = def.max_node_energy;
    interpol_def.romberg = def.order_of_interpolation;
    interpol_def.rational = true;
    interpol_def.isLog = true;
    interpol_def.rombergY = def.order_of_interpolation;
    interpol_def.logSubst = true;

    unique_ptr<InterpolantBuilder> builder(
        new Interpolant1DBuilder(interpol_def));

    auto aux = Helper::InitializeInterpolation(
        "dEdx", std::move(builder), GetHash(), def);

    return aux;
}

unique_ptr<Interpolant> CrossSectionInterpolant::init_de2dx_interpolation(
    const InterpolationDef& def)
{
    Interpolant1DBuilder::Definition interpol_def;
    interpol_def.function1d = [this](double energy) {
        return CrossSectionIntegral::CalculatedE2dx(energy);
    };
    interpol_def.max = def.nodes_continous_randomization;
    interpol_def.xmin = parametrization_->GetLowerEnergyLim();
    interpol_def.xmax = def.max_node_energy;
    interpol_def.romberg = def.order_of_interpolation;
    interpol_def.isLog = true;
    interpol_def.rombergY = def.order_of_interpolation;

    unique_ptr<InterpolantBuilder> builder(
        new Interpolant1DBuilder(interpol_def));

    return Helper::InitializeInterpolation(
        "dE2dx", std::move(builder), GetHash(), def);
}

vector<unique_ptr<Interpolant>>
CrossSectionInterpolant::init_dndx_interpolation(const InterpolationDef& def)
{
    Interpolant2DBuilder::Definition interpol_def;
    interpol_def.max1 = def.nodes_cross_section;
    interpol_def.x1min = parametrization_->GetLowerEnergyLim();
    interpol_def.x1max = def.max_node_energy;
    interpol_def.max2 = def.nodes_cross_section;
    interpol_def.x2min = 0.0;
    interpol_def.x2max = 1.0;
    interpol_def.romberg1 = def.order_of_interpolation;
    interpol_def.isLog1 = true;
    interpol_def.romberg2 = def.order_of_interpolation;
    interpol_def.rombergY = def.order_of_interpolation;
    interpol_def.rationalY = true;

    vector<unique_ptr<InterpolantBuilder>> builder;
    for (const auto& dndx : dndx_integral_) {
        interpol_def.function2d = dndx;
        builder.emplace_back(new Interpolant2DBuilder(interpol_def));
    }

    return Helper::InitializeInterpolation(
        "dNdx", builder, GetHash(), def);
}

double CrossSectionInterpolant::CalculatedEdx(double energy)
{
    return dedx_interpolant_->Interpolate(energy);
}

double CrossSectionInterpolant::CalculatedE2dx(double energy)
{
    return de2dx_interpolant_->Interpolate(energy);
}

vector<double> CrossSectionInterpolant::CalculatedNdx(double energy)
{
    vector<double> rates;
    for (auto& interpol : dndx_interpolants_)
        rates.emplace_back(interpol->Interpolate(energy, 1.));
    return rates;
}

vector<double> CrossSectionInterpolant::CalculateStochasticLoss(
    double energy, const vector<double>& component_rates)
{
    vector<double> relativ_losses;
    auto sampled_component_rate = component_rates.cbegin();
    for (auto& interpol : dndx_interpolants_)
        relativ_losses.emplace_back(interpol->FindLimit(energy, *sampled_component_rate));

    auto relativ_loss = relativ_losses.begin();
    for (const auto& comp : parametrization_->GetComponents()) {
        parametrization_->SetCurrentComponent(comp);
        *relativ_loss = transform_relativ_loss(*relativ_loss, energy);
        ++relativ_loss;
    }

    return relativ_losses;
}
