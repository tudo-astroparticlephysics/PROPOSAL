
#include <cmath>
#include <functional>

#include "PROPOSAL/crossection/CrossSectionInterpolant.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/medium/Medium.h"

using namespace PROPOSAL;
using std::bind;
using std::placeholders::_1;
using std::placeholders::_2;

double CrossSectionInterpolant::logarithm_trafo(
    double v, double v_cut, double v_max) const
{
    return v_cut * std::exp(v * std::log(v_max / v_cut));
}

unique_ptr<Interpolant> CrossSectionInterpolant::init_dedx_interpolation(
    const InterpolationDef& def)
{
    Interpolant1DBuilder::Definition interpol_def;
    interpol_def.function1d
        = bind(&CrossSectionIntegral::CalculatedEdx, this, _1);
    interpol_def.max = def.nodes_cross_section;
    interpol_def.xmin = parametrization_->GetLowerEnergyLim();
    interpol_def.xmax = def.max_node_energy;
    interpol_def.romberg = def.order_of_interpolation;
    interpol_def.rational = true;
    interpol_def.isLog = true;
    interpol_def.rombergY = def.order_of_interpolation;
    interpol_def.logSubst = true;

    Interpolant1DBuilder builder(interpol_def);
    auto hash = parametrization_->GetHash();

    return Helper::InitializeInterpolation("dEdx", builder, hash, def);
}

unique_ptr<Interpolant> CrossSectionInterpolant::init_de2dx_interpolation(
    const InterpolationDef& def)
{
    Interpolant1DBuilder::Definition interpol_def;
    interpol_def.function1d
        = bind(&CrossSectionIntegral::CalculatedE2dx, this, _1);
    interpol_def.max = def.nodes_continous_randomization;
    interpol_def.xmin = parametrization_->GetLowerEnergyLim();
    interpol_def.xmax = def.max_node_energy;
    interpol_def.romberg = def.order_of_interpolation;
    interpol_def.isLog = true;
    interpol_def.rombergY = def.order_of_interpolation;

    Interpolant1DBuilder builder(interpol_def);
    auto hash = parametrization_->GetHash();

    return Helper::InitializeInterpolation("dE2dx", builder, hash, def);
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
        "dNdx", builder, parametrization_->GetHash(), def);
}

double CrossSectionInterpolant::CalculatedEdx(double energy)
{
    return dedx_interpolant_->Interpolate(energy);
}

double CrossSectionInterpolant::CalculatedE2dx(double energy)
{
    return de2dx_interpolant_->Interpolate(energy);
}

double CrossSectionInterpolant::CalculatedNdx(double energy, double rnd)
{
    vector<double> rates;
    for (auto& interpol : dndx_interpolants_)
        rates.emplace_back(interpol->Interpolate(energy, 1.));

    auto total_rate = accumulate(rates.begin(), rates.end(), 0.);

    size_t nth_component = 0;
    for (const auto& rate : rates) {
        total_rate -= rate / rnd;
        if (total_rate < 0)
            return dndx_interpolants_.at(nth_component)
                ->Interpolate(energy, total_rate);
        nth_component += 1;
    }

    return total_rate;
}

size_t CrossSectionInterpolant::GetHash() const
{
    size_t hash_digest = 0;
    hash_combine(hash_digest, hash_interpol_def, CrossSection::GetHash());

    return hash_digest;
}
