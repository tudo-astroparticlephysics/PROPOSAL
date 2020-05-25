
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

namespace PROPOSAL {

template <class Param> class CrossSectionInterpolant : public CrossSection {
    InterpolationDef def;
    Param param;

    using base_param_t = typename std::add_lvalue_reference<
        typename decay<Param>::type::base_param_t>::type;

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
        const ParticleDef&, const Component&, double, double);

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
    return Helper::InitializeInterpolation(
        "dE2dx", std::move(builder), hash, def);
}

template <typename Param>
unique_ptr<Interpolant> build_dndx(Param&& param, const ParticleDef& p_def,
    const Component& comp, const EnergyCutSettings& cut,
    const InterpolationDef& def, size_t hash)
{
    Integral integral;
    Interpolant2DBuilder::Definition interpol_def;
    interpol_def.max1 = def.nodes_cross_section;
    interpol_def.x1min = param.GetLowerEnergyLim(p_def);
    std::cout << "interpol_def.x1min: " << interpol_def.x1min << std::endl;
    interpol_def.x1max = def.max_node_energy;
    interpol_def.max2 = def.nodes_cross_section;
    interpol_def.x2min = 0.0;
    interpol_def.x2max = 1.0;
    interpol_def.romberg1 = def.order_of_interpolation;
    interpol_def.isLog1 = true;
    interpol_def.romberg2 = def.order_of_interpolation;
    interpol_def.rombergY = def.order_of_interpolation;
    interpol_def.rationalY = true;
    interpol_def.function2d = [&integral, &param, &p_def, &comp, &cut](
                                  double energy, double v) {
        auto physical_lim = param.GetKinematicLimits(p_def, comp, energy);
        auto v_cut = cut.GetCut(physical_lim, energy);
        v = transform_relativ_loss(v_cut, physical_lim.vMax, v);
        return integrate_dndx(integral, param, p_def, comp, energy, v_cut, v);
    };
    auto builder = make_unique<Interpolant2DBuilder>(interpol_def);
    return Helper::InitializeInterpolation(
        "dNdx", std::move(builder), hash, def);
}

template <class Param>
double CrossSectionInterpolant<Param>::CalculatedEdx(
    const ParticleDef& p_def, const Medium& medium, double energy)
{
    auto hash = GetHash(p_def, medium);
    auto search = dedx_interpolants.find(hash);
    if (search != dedx_interpolants.end())
        return search->second->Interpolate(energy);
    dedx_interpolants[hash] = build_dedx(reinterpret_cast<base_param_t>(param),
        p_def, medium, *cuts_, def, hash);
    return dedx_interpolants[hash]->Interpolate(energy);
}

template <class Param>
double CrossSectionInterpolant<Param>::CalculatedE2dx(
    const ParticleDef& p_def, const Medium& medium, double energy)
{
    auto hash = GetHash(p_def, medium);
    auto search = de2dx_interpolants.find(hash);
    if (search != de2dx_interpolants.end())
        return search->second->Interpolate(energy);
    de2dx_interpolants[hash]
        = build_de2dx(reinterpret_cast<base_param_t>(param), p_def, medium,
            *cuts_, def, hash);
    return de2dx_interpolants[hash]->Interpolate(energy);
}

template <class Param>
unordered_map<size_t, double> CrossSectionInterpolant<Param>::CalculatedNdx(
    const ParticleDef& p_def, const Medium& medium, double energy)
{
    unordered_map<size_t, double> rates;
    for (const auto comp : medium.GetComponents()) {
        auto hash = GetHash(p_def, comp);
        auto search = dndx_interpolants.find(hash);
        if (search != dndx_interpolants.end())
            rates[hash] = search->second->Interpolate(energy, 1.);
        else {
            dndx_interpolants[hash]
                = build_dndx(reinterpret_cast<base_param_t>(param), p_def, comp,
                    *cuts_, def, hash);
            rates[hash] = dndx_interpolants[hash]->Interpolate(energy, 1.);
        }
    }
    return rates;
}

template <class Param>
double CrossSectionInterpolant<Param>::CalculateStochasticLoss(
    const ParticleDef& p_def, const Component& comp, double energy, double rate)
{
    auto hash = GetHash(p_def, comp);
    auto v = dndx_interpolants[hash]->FindLimit(energy, rate);
    auto physical_lim = param.GetKinematicLimits(p_def, comp, energy);
    auto v_cut = cuts_->GetCut(physical_lim, energy);
    return transform_relativ_loss(v_cut, physical_lim.vMax, v);
}

template <typename Param>
size_t CrossSectionInterpolant<Param>::GetHash(
    const ParticleDef& p_def, const Medium& medium) const
{
    size_t hash_digest = 0;
    hash_combine(hash_digest, param.GetHash(), p_def.GetHash(),
        medium.GetHash(), def.GetHash());
    if (cuts_)
        hash_combine(hash_digest, cuts_->GetHash());
    return hash_digest;
}

template <typename Param>
size_t CrossSectionInterpolant<Param>::GetHash(
    const ParticleDef& p_def, const Component& comp) const
{
    size_t hash_digest = 0;
    hash_combine(hash_digest, param.GetHash(), p_def.GetHash(), comp.GetHash(),
        def.GetHash());
    if (cuts_)
        hash_combine(hash_digest, cuts_->GetHash());
    return hash_digest;
}

template <typename Param>
double CrossSectionInterpolant<Param>::GetLowerEnergyLim(
    const ParticleDef& p_def) const
{
    return param.GetLowerEnergyLim(p_def);
}

} // namespace PROPOSAL
