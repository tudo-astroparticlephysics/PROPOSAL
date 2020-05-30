
/******************************************************************************
 *                                                                            *
 * This file is part of the simulation tool PROPOSAL.                         *
 *                                                                            *
 * Copyright (C) 2017 TU Dortmund University, Department of Physics,          *
 *                    Chair Experimental Physics 5b                           *
 *                                                                            *
 * This software may be modified and distributed under the terms of a         *
 * modified GNU Lesser General Public Licence version 3 (LGPL),               *
 * copied verbatim in the file "LICENSE".                                     *
 *                                                                            *
 * Modifcations to the LGPL License:                                          *
 *                                                                            *
 *      1. The user shall acknowledge the use of PROPOSAL by citing the       *
 *         following reference:                                               *
 *                                                                            *
 *         J.H. Koehne et al.  Comput.Phys.Commun. 184 (2013) 2070-2090 DOI:  *
 *         10.1016/j.cpc.2013.04.001                                          *
 *                                                                            *
 *      2. The user should report any bugs/errors or improvments to the       *
 *         current maintainer of PROPOSAL or open an issue on the             *
 *         GitHub webpage                                                     *
 *                                                                            *
 *         "https://github.com/tudo-astroparticlephysics/PROPOSAL"            *
 *                                                                            *
 ******************************************************************************/

#pragma once

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/CrossSectionIntegral.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/methods.h"

#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/particle/ParticleDef.h"

using std::get;

namespace PROPOSAL {

template <class Param> class CrossSectionInterpolant : public CrossSection {
    InterpolationDef def;
    Param param;

    using base_param_t = typename decay<Param>::type::base_param_t;
    using base_param_ref_t = typename add_lvalue_reference<base_param_t>::type;

    rates_t CalculatedNdx_impl(const ParticleDef& p_def, const Medium& medium,
        double energy, std::true_type);
    rates_t CalculatedNdx_impl(const ParticleDef& p_def, const Medium& medium,
        double energy, std::false_type);

    double CalculateStochasticLoss_impl(const ParticleDef&, const Medium&,
        const Component&, double, double, std::false_type);
    double CalculateStochasticLoss_impl(const ParticleDef&, const Medium&,
        const Component&, double, double, std::true_type);

protected:
    unordered_map<size_t, unique_ptr<Interpolant>> dedx_interpolants;
    unordered_map<size_t, unique_ptr<Interpolant>> de2dx_interpolants;
    unordered_map<size_t, unique_ptr<Interpolant>> dndx_interpolants;

public:
    CrossSectionInterpolant(
        Param&&, shared_ptr<const EnergyCutSettings>, const InterpolationDef&);

    double CalculatedEdx(const ParticleDef&, const Medium&, double) override;
    double CalculatedE2dx(const ParticleDef&, const Medium&, double) override;
    rates_t CalculatedNdx(const ParticleDef&, const Medium&, double) override;
    double CalculateStochasticLoss(
        const ParticleDef&, const Medium&, const Component&, double, double);

    double GetLowerEnergyLim(const ParticleDef&) const override;
    size_t GetHash(const ParticleDef&, const Medium&) const override;
    size_t GetHash(const ParticleDef&, const Component&) const override;
};

double transform_relativ_loss(double v_cut, double v_max, double v);

template <class Param>
CrossSectionInterpolant<Param>::CrossSectionInterpolant(Param&& param,
    shared_ptr<const EnergyCutSettings> cut, const InterpolationDef& def)
    : CrossSection(cut)
    , param(param)
    , def(def)
{
}

template <typename Param>
unique_ptr<Interpolant> build_dedx(Param&& param, const ParticleDef& p_def,
    const Medium& medium, const EnergyCutSettings& cut,
    const InterpolationDef& def, size_t hash)
{
    Integral integral;
    Interpolant1DBuilder::Definition interpol_def;
    interpol_def.function1d
        = [&integral, &param, &p_def, &medium, &cut](double energy) {
              return calculate_dedx(param, integral, p_def, medium, cut, energy,
                  typename decay<Param>::type::only_stochastic{});
          };
    interpol_def.max = def.nodes_cross_section;
    interpol_def.xmin = param.GetLowerEnergyLim(p_def);
    interpol_def.xmax = def.max_node_energy;
    interpol_def.romberg = def.order_of_interpolation;
    interpol_def.rational = true;
    interpol_def.isLog = true;
    interpol_def.rombergY = def.order_of_interpolation;
    interpol_def.logSubst = true;
    hash_combine(hash, param.GetHash(), cut.GetHash());
    return Helper::InitializeInterpolation(
        "dEdx", make_unique<Interpolant1DBuilder>(interpol_def), hash, def);
}

template <typename Param>
unique_ptr<Interpolant> build_de2dx(Param&& param, const ParticleDef& p_def,
    const Medium& medium, const EnergyCutSettings& cut,
    const InterpolationDef& def, size_t hash)
{
    Integral integral;
    Interpolant1DBuilder::Definition interpol_def;
    interpol_def.function1d
        = [&integral, &param, &p_def, &medium, &cut](double energy) {
              return calculate_de2dx(param, integral, p_def, medium, cut,
                  energy, typename decay<Param>::type::only_stochastic{});
          };
    interpol_def.max = def.nodes_continous_randomization;
    interpol_def.xmin = param.GetLowerEnergyLim(p_def);
    interpol_def.xmax = def.max_node_energy;
    interpol_def.romberg = def.order_of_interpolation;
    interpol_def.isLog = true;
    interpol_def.rombergY = def.order_of_interpolation;
    auto builder = make_unique<Interpolant1DBuilder>(interpol_def);
    hash_combine(hash, param.GetHash(), cut.GetHash());
    return Helper::InitializeInterpolation(
        "dE2dx", std::move(builder), hash, def);
}

template <typename P, typename M>
unique_ptr<Interpolant> build_dndx(P&& param, const ParticleDef& p_def,
     M&& medium, const EnergyCutSettings& cut,
    const InterpolationDef& def, size_t hash)
{
    Integral integral;
    Interpolant2DBuilder::Definition interpol_def;
    interpol_def.max1 = def.nodes_cross_section;
    interpol_def.x1min = param.GetLowerEnergyLim(p_def);
    interpol_def.x1max = def.max_node_energy;
    interpol_def.max2 = def.nodes_cross_section;
    interpol_def.x2min = 0.0;
    interpol_def.x2max = 1.0;
    interpol_def.romberg1 = def.order_of_interpolation;
    interpol_def.isLog1 = true;
    interpol_def.romberg2 = def.order_of_interpolation;
    interpol_def.rombergY = def.order_of_interpolation;
    interpol_def.rationalY = true;
    interpol_def.function2d = [&integral, &param, &p_def, &medium, &cut](
                                  double energy, double v) {
        auto lim = param.GetKinematicLimits(p_def, medium, energy);
        auto v_cut = cut.GetCut(lim, energy);
        v = transform_relativ_loss(v_cut, get<Parametrization::V_MAX>(lim), v);
        return integrate_dndx(integral, param, p_def, medium, energy, v_cut, v);
    };
    auto builder = make_unique<Interpolant2DBuilder>(interpol_def);
    hash_combine(hash, param.GetHash(), cut.GetHash());
    return Helper::InitializeInterpolation(
        "dNdx", std::move(builder), hash, def);
}

template <class Param>
double CrossSectionInterpolant<Param>::CalculatedEdx(
    const ParticleDef& p_def, const Medium& medium, double energy)
{
    auto hash = size_t{ 0 };
    hash_combine(hash, p_def.GetHash(), medium.GetHash());
    auto search = dedx_interpolants.find(hash);
    if (search != dedx_interpolants.end())
        return search->second->Interpolate(energy);
    dedx_interpolants[hash]
        = build_dedx(reinterpret_cast<base_param_ref_t>(param), p_def, medium,
            *cuts_, def, hash);
    return dedx_interpolants[hash]->Interpolate(energy);
}

template <class Param>
double CrossSectionInterpolant<Param>::CalculatedE2dx(
    const ParticleDef& p_def, const Medium& medium, double energy)
{
    auto hash = size_t{ 0 };
    hash_combine(hash, p_def.GetHash(), medium.GetHash());
    auto search = de2dx_interpolants.find(hash);
    if (search != de2dx_interpolants.end())
        return search->second->Interpolate(energy);
    de2dx_interpolants[hash]
        = build_de2dx(reinterpret_cast<base_param_ref_t>(param), p_def, medium,
            *cuts_, def, hash);
    return de2dx_interpolants[hash]->Interpolate(energy);
}

template <typename Param>
rates_t CrossSectionInterpolant<Param>::CalculatedNdx_impl(
    const ParticleDef& p_def, const Medium& medium, double energy,
    std::true_type)
{
    rates_t rates;
    for (auto& c : medium.GetComponents()) {
        auto hash = size_t{ 0 };
        hash_combine(hash, p_def.GetHash(), c.GetHash());
        auto search = dndx_interpolants.find(hash);
        if (search != dndx_interpolants.end())
            rates[&c] = search->second->Interpolate(energy, 1.);
        else {
            dndx_interpolants[hash]
                = build_dndx(reinterpret_cast<base_param_ref_t>(param), p_def,
                    c, *cuts_, def, hash);
            rates[&c] = dndx_interpolants[hash]->Interpolate(energy, 1.);
        }
    }
    return rates;
}

template <typename Param>
rates_t CrossSectionInterpolant<Param>::CalculatedNdx_impl(
    const ParticleDef& p_def, const Medium& medium, double energy,
    std::false_type)
{
    rates_t rates;
    auto hash = size_t{ 0 };
    hash_combine(hash, p_def.GetHash(), medium.GetHash());
    auto search = dndx_interpolants.find(hash);
    if (search != dndx_interpolants.end())
        rates[nullptr] = search->second->Interpolate(energy, 1.);
    else {
        dndx_interpolants[hash]
            = build_dndx(reinterpret_cast<base_param_ref_t>(param), p_def,
                medium, *cuts_, def, hash);
        rates[nullptr] = dndx_interpolants[hash]->Interpolate(energy, 1.);
    }
    return rates;
}

template <class Param>
rates_t CrossSectionInterpolant<Param>::CalculatedNdx(
    const ParticleDef& p_def, const Medium& medium, double energy)
{
    return CalculatedNdx_impl(
        p_def, medium, energy, typename base_param_t::component_wise{});
}

template <typename Param>
double CrossSectionInterpolant<Param>::CalculateStochasticLoss_impl(
    const ParticleDef& p_def, const Medium& medium, const Component& comp,
    double energy, double rate, std::true_type)
{
    auto hash = size_t{ 0 };
    hash_combine(hash, p_def.GetHash(), comp.GetHash());
    auto v = dndx_interpolants[hash]->FindLimit(energy, rate);
    auto lim = param.GetKinematicLimits(p_def, comp, energy);
    auto v_cut = cuts_->GetCut(lim, energy);
    return transform_relativ_loss(v_cut, get<Parametrization::V_MAX>(lim), v);
}

template <typename Param>
double CrossSectionInterpolant<Param>::CalculateStochasticLoss_impl(
    const ParticleDef& p_def, const Medium& medium, const Component& comp,
    double energy, double rate, std::false_type)
{
    auto hash = size_t{ 0 };
    hash_combine(hash, p_def.GetHash());
    auto v = dndx_interpolants[hash]->FindLimit(energy, rate);
    auto lim = param.GetKinematicLimits(p_def, medium, energy);
    auto v_cut = cuts_->GetCut(lim, energy);
    return transform_relativ_loss(v_cut, get<Parametrization::V_MAX>(lim), v);
}

template <class Param>
double CrossSectionInterpolant<Param>::CalculateStochasticLoss(
    const ParticleDef& p_def, const Medium& medium, const Component& comp,
    double energy, double rate)
{
    return CalculateStochasticLoss_impl(p_def, medium, comp, energy, rate,
        typename base_param_t::component_wise{});
}

template <typename Param>
size_t CrossSectionInterpolant<Param>::GetHash(
    const ParticleDef& p_def, const Medium& medium) const
{
    auto hash = size_t{ 0 };
    hash_combine(hash, p_def.GetHash(), medium.GetHash(), def.GetHash());
    if (cuts_)
        hash_combine(hash, cuts_->GetHash());
    return hash;
}

template <typename Param>
size_t CrossSectionInterpolant<Param>::GetHash(
    const ParticleDef& p_def, const Component& comp) const
{
    auto hash = size_t{ 0 };
    hash_combine(hash, p_def.GetHash(), comp.GetHash(), def.GetHash());
    if (cuts_)
        hash_combine(hash, cuts_->GetHash());
    return hash;
}

template <typename Param>
double CrossSectionInterpolant<Param>::GetLowerEnergyLim(
    const ParticleDef& p_def) const
{
    return param.GetLowerEnergyLim(p_def);
}

} // namespace PROPOSAL
