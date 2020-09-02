
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
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXBuilder.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/methods.h"

#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/particle/ParticleDef.h"

using std::get;
using std::unordered_map;

namespace PROPOSAL {

template <typename Param, typename P, typename M>
class CrossSectionInterpolant : public crosssection_t<P, M> {

    using param_t = typename decay<Param>::type;
    using particle_t = typename decay<P>::type;
    using medium_t = typename decay<M>::type;
    using base_param_ref_t =
        typename add_lvalue_reference<typename param_t::base_param_t>::type;

    param_t param;
    particle_t p_def;
    medium_t medium;
    shared_ptr<const EnergyCutSettings> cut;
    InterpolationDef def;

    dndx_map_t dndx_map;

    double CalculateStochasticLoss_impl(
        const Component&, double, double, std::true_type, std::false_type);
    double CalculateStochasticLoss_impl(
        const Component&, double, double, std::false_type, std::false_type);
    double CalculateStochasticLoss_impl(
            const Component&, double, double, bool, std::true_type);
protected:
    unique_ptr<Interpolant> dedx;
    unique_ptr<Interpolant> de2dx;
    /* unordered_map<const Component*, unique_ptr<Interpolant>> dndx; */

public:
    CrossSectionInterpolant(Param&& _param, P&& _p_def, M&& _medium,
        shared_ptr<const EnergyCutSettings> _cut, const InterpolationDef& _def)
        : CrossSection<typename decay<P>::type, typename decay<M>::type>()
        , param(std::forward<Param>(_param)) // only for back transformation
        , p_def(std::forward<P>(_p_def))     // needed TODO: Maximilian Sackel
        , medium(std::forward<M>(_medium))   // 2 Jun. 2020
        , cut(_cut)
        , def(_def)
        , dndx_map(build_cross_section_dndx(param, p_def, medium, cut, true))
    {
        if (typename param_t::only_stochastic{} == true and cut != nullptr) {
            throw std::invalid_argument("CrossSections of parametrizations that are only stochastic do not use"
                                        "EnergyCuts. Pass a nullptr as an EnergyCut instead.");
        }

        if (cut != nullptr) {
            // Only for a defined EnergyCut, dEdx and dE2dx return non-zero values
            dedx = build_dedx(reinterpret_cast<base_param_ref_t>(param), p_def,
                              medium, *cut, def, 0);
            de2dx = build_de2dx(reinterpret_cast<base_param_ref_t>(param), p_def,
                                medium, *cut, def, 0);
        }
    }

    inline double CalculatedEdx(double energy) override
    {
        if (dedx == nullptr)
            return 0;
        return dedx->Interpolate(energy);
    }

    inline double CalculatedE2dx(double energy) override
    {
        if (de2dx == nullptr)
            return 0;
        return de2dx->Interpolate(energy);
    }
    inline rates_t CalculatedNdx(double energy) override
    {
        auto rates = rates_t();
        for (auto& dndx : dndx_map) {
            //TODO: dNdx interpolant results for individual components can become negative for small energies
            // Instead of clipping these values to zero, the interpolant should be revised (jm)
            rates[dndx.first] = std::max(dndx.second->Calculate(energy), 0.);
            if (dndx.first)
                rates[dndx.first] /= medium.GetSumNucleons()
                    / (dndx.first->GetAtomInMolecule()
                          * dndx.first->GetAtomicNum());
        }
        return rates;
    }
    inline double CalculateStochasticLoss(
        const Component& comp, double energy, double rate) override
    {
        return CalculateStochasticLoss_impl(comp, energy, rate,
            typename param_t::base_param_t::component_wise{},
            typename param_t::base_param_t::only_stochastic{});
    }
    inline size_t GetHash() const noexcept override
    {
        auto hash_digest = size_t{ 0 };
        hash_combine(hash_digest, param.GetHash(), p_def.GetHash(),
            medium.GetHash(), def.GetHash());
        if (cut != nullptr)
            hash_combine(hash_digest, cut->GetHash());
        return hash_digest;
    }
    inline double GetLowerEnergyLim() const override
    {
        return param.GetLowerEnergyLim(p_def);
    }
    inline InteractionType GetInteractionType() const noexcept override
    {
        return param.interaction_type;
    }
};

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
                  typename decay<Param>::type::component_wise{});
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
                  energy, typename decay<Param>::type::component_wise{});
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

/* using dndx_func_t = function<double(double, double)>; */

/* template <typename Param> */
/* unordered_map<const Component*, dndx_func_t> build_dndx_functions(Param&&
 * param, */
/*     const ParticleDef& p_def, const Medium& medium, */
/*     const EnergyCutSettings& cut, std::true_type) */
/* { */
/*     unordered_map<const Component*, dndx_func_t> dndx_functions; */
/*     for (const auto& comp : medium.GetComponents()) { */
/*         dndx_functions[&comp] */
/*             = [&param, &p_def, &comp, &cut](double energy, double v) { */
/*                   Integral integral; */
/*                   auto lim = param.GetKinematicLimits(p_def, comp, energy);
 */
/*                   auto v_cut = cut.GetCut(lim, energy); */
/*                   v = transform_relativ_loss( */
/*                       v_cut, get<Parametrization::V_MAX>(lim), v); */
/*                   return integrate_dndx( */
/*                       integral, param, p_def, comp, energy, v_cut, v); */
/*               }; */
/*     } */
/*     return dndx_functions; */
/* } */

/* template <typename Param> */
/* unordered_map<const Component*, dndx_func_t> build_dndx_functions(Param&&
 * param, */
/*     const ParticleDef& p_def, const Medium& medium, */
/*     const EnergyCutSettings& cut, std::false_type) */
/* { */
/*     Integral integral; */
/*     unordered_map<const Component*, dndx_func_t> dndx_functions; */
/*     auto dndx_func = [&param, &p_def, &medium, &cut](double energy, double v)
 * { */
/*         Integral integral; */
/*         auto lim = param.GetKinematicLimits(p_def, medium, energy); */
/*         auto v_cut = cut.GetCut(lim, energy); */
/*         v = transform_relativ_loss(v_cut, get<Parametrization::V_MAX>(lim),
 * v); */
/*         return integrate_dndx(integral, param, p_def, medium, energy, v_cut,
 * v); */
/*     }; */
/*     dndx_functions[nullptr] = dndx_func; */
/*     return dndx_functions; */
/* } */

/* template <typename Param> */
/* unordered_map<const Component*, unique_ptr<Interpolant>> build_dndx( */
/*     Param&& param, const ParticleDef& p_def, const Medium& medium, */
/*     const EnergyCutSettings& cut, const InterpolationDef& def, size_t hash)
 */
/* { */
/*     Interpolant2DBuilder::Definition interpol_def; */
/*     interpol_def.max1 = def.nodes_cross_section; */
/*     interpol_def.x1min = param.GetLowerEnergyLim(p_def); */
/*     interpol_def.x1max = def.max_node_energy; */
/*     interpol_def.max2 = def.nodes_cross_section; */
/*     interpol_def.x2min = 0.0; */
/*     interpol_def.x2max = 1.0; */
/*     interpol_def.romberg1 = def.order_of_interpolation; */
/*     interpol_def.isLog1 = true; */
/*     interpol_def.romberg2 = def.order_of_interpolation; */
/*     interpol_def.rombergY = def.order_of_interpolation; */
/*     interpol_def.rationalY = true; */

/*     using param_t = typename decay<Param>::type; */
/*     unordered_map<const Component*, dndx_func_t> dndx_functions */
/*         = build_dndx_functions( */
/*             param, p_def, medium, cut, typename param_t::component_wise{});
 */

/*     unordered_map<const Component*, unique_ptr<Interpolant>> dndx_interpol;
 */
/*     for (const auto& dndx_func : dndx_functions) { */
/*         interpol_def.function2d = dndx_func.second; */
/*         auto builder = make_unique<Interpolant2DBuilder>(interpol_def); */
/*         hash_combine(hash, param.GetHash(), cut.GetHash()); */
/*         dndx_interpol[dndx_func.first] = Helper::InitializeInterpolation( */
/*             "dNdx", std::move(builder), hash, def); */
/*     } */
/*     return dndx_interpol; */
/* } */

template <typename Param, typename P, typename M>
double CrossSectionInterpolant<Param, P, M>::CalculateStochasticLoss_impl(
    const Component& comp, double energy, double rate, std::true_type, std::false_type)
{
    auto weight_for_rate_in_medium = medium.GetSumNucleons()
        / (comp.GetAtomInMolecule() * comp.GetAtomicNum());
    return dndx_map[&comp]->GetUpperLimit(
        energy, rate * weight_for_rate_in_medium);
}

template <typename Param, typename P, typename M>
double CrossSectionInterpolant<Param, P, M>::CalculateStochasticLoss_impl(
    const Component&, double energy, double rate, std::false_type, std::false_type)
{
    return dndx_map[nullptr]->GetUpperLimit(energy, rate);
}

template <typename Param, typename P, typename M>
double CrossSectionInterpolant<Param, P, M>::CalculateStochasticLoss_impl(
    const Component&, double energy, double, bool, std::true_type)
{
    return energy;
}

} // namespace PROPOSAL
