
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
    , v_trafo(init_v_trafo())
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

    auto builder = make_unique<Interpolant1DBuilder>(interpol_def);

    return Helper::InitializeInterpolation(
        "dEdx", std::move(builder), GetHash(), def);
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

    auto builder = make_unique<Interpolant1DBuilder>(interpol_def);

    return Helper::InitializeInterpolation(
        "dE2dx", std::move(builder), GetHash(), def);
}

unordered_map<size_t, unique_ptr<Interpolant>>
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

    unordered_map<size_t, unique_ptr<Interpolant>> dndx_;
    for (const auto& comp : parametrization_->GetComponents()) {
        interpol_def.function2d = dndx_integral_[comp.GetHash()];

        auto builder = make_unique<Interpolant2DBuilder>(interpol_def);
        auto interpolant = Helper::InitializeInterpolation(
            "dNdx", std::move(builder), GetHash(), def);

        dndx_[comp.GetHash()] = std::move(interpolant);
    }

    return dndx_;
}

unordered_map<size_t, function<double(double, double)>>
CrossSectionInterpolant::init_v_trafo()
{
    unordered_map<size_t, function<double(double, double)>> trafo;
    for (const auto& comp : parametrization_->GetComponents()) {
        parametrization_->SetCurrentComponent(comp);
        trafo[comp.GetHash()] = bind(
            &CrossSectionInterpolant::transform_relativ_loss, this, _1, _2);
    }
    return trafo;
}

double CrossSectionInterpolant::CalculatedEdx(double energy)
{
    return dedx_interpolant_->Interpolate(energy);
}

double CrossSectionInterpolant::CalculatedE2dx(double energy)
{
    return de2dx_interpolant_->Interpolate(energy);
}

unordered_map<size_t, double> CrossSectionInterpolant::CalculatedNdx(
    double energy)
{
    unordered_map<size_t, double> rates;
    for (auto& comp : parametrization_->GetComponents())
        rates[comp.GetHash()]
            = dndx_interpolants_[comp.GetHash()]->Interpolate(energy, 1.);

    return rates;
}

double CrossSectionInterpolant::CalculateStochasticLoss(
    double energy, double rate, size_t comp_hash)
{
    auto v = dndx_interpolants_[comp_hash]->FindLimit(energy, rate);
    return v_trafo[comp_hash](v, energy);
}
