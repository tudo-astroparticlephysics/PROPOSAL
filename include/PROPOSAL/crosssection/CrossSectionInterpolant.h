
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
#include "PROPOSAL/crosssection/CrossSectionInterpolantBase.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXBuilder.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/methods.h"

#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/particle/ParticleDef.h"


namespace PROPOSAL {

template <typename Param, typename P, typename M>
class CrossSectionInterpolant : public crosssection_t<P, M>, public CrossSectionInterpolantBase {

    using param_t = typename std::decay<Param>::type;
    using particle_t = typename std::decay<P>::type;
    using medium_t = typename std::decay<M>::type;
    using base_param_ref_t =
        typename std::add_lvalue_reference<typename param_t::base_param_t>::type;

    param_t param;
    particle_t p_def;
    medium_t medium;
    shared_ptr<const EnergyCutSettings> cut;


    dndx_map_t dndx_map;

    double CalculateStochasticLoss_impl(
        std::shared_ptr<const Component> const&, double, double, std::true_type, std::false_type);
    double CalculateStochasticLoss_impl(
        std::shared_ptr<const Component> const&, double, double, std::false_type, std::false_type);
    double CalculateStochasticLoss_impl(
            std::shared_ptr<const Component> const&, double, double, bool, std::true_type);
protected:
    unique_ptr<Interpolant> dedx;
    unique_ptr<Interpolant> de2dx;
    /* unordered_map<const Component*, unique_ptr<Interpolant>> dndx; */

public:
    CrossSectionInterpolant(Param&& _param, P&& _p_def, M&& _medium,
        shared_ptr<const EnergyCutSettings> _cut)
        : CrossSection<typename std::decay<P>::type, typename std::decay<M>::type>()
        , param(std::forward<Param>(_param)) // only for back transformation
        , p_def(std::forward<P>(_p_def))     // needed TODO: Maximilian Sackel
        , medium(std::forward<M>(_medium))   // 2 Jun. 2020
        , cut(_cut)
        , dndx_map(build_cross_section_dndx(param, p_def, medium, cut, true))
    {
        if (typename param_t::only_stochastic{} == true and cut != nullptr) {
            throw std::invalid_argument("CrossSections of parametrizations that are only stochastic do not use"
                                        "EnergyCuts. Pass a nullptr as an EnergyCut instead.");
        }

        if (cut != nullptr) {
            // Only for a defined EnergyCut, dEdx and dE2dx return non-zero values
            dedx = build_dedx(reinterpret_cast<base_param_ref_t>(param), p_def,
                              medium, *cut, dEdx_def);
            de2dx = build_de2dx(reinterpret_cast<base_param_ref_t>(param), p_def,
                                medium, *cut, dE2dx_def);
        }
    }
    std::vector<std::shared_ptr<const Component>> GetTargets() const noexcept final
    {
        std::vector<std::shared_ptr<const Component>> targets;
        for (auto& dndx : dndx_map)
        {
            targets.emplace_back(dndx.first);
        }
        return targets;
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
    // inline rates_t CalculatedNdx(double energy) override
    // {
    //     auto rates = rates_t();
    //     for (auto& dndx : dndx_map) {
    //         //TODO: dNdx interpolant results for individual components can become negative for small energies
    //         // Instead of clipping these values to zero, the interpolant should be revised (jm)
    //         rates[dndx.first] = std::max(dndx.second->Calculate(energy), 0.);
    //         if (dndx.first)
    //             rates[dndx.first] /= medium.GetSumNucleons()
    //                 / (dndx.first->GetAtomInMolecule()
    //                       * dndx.first->GetAtomicNum());
    //     }
    //     return rates;
    // }
    inline double CalculatedNdx(double energy, std::shared_ptr<const Component> comp_ptr) override
    {
        //TODO: dNdx interpolant results for individual components can become negative for small energies
        // Instead of clipping these values to zero, the interpolant should be revised (jm)
        if (comp_ptr == nullptr)
        {
            double dndx_all = 0.;
            double tmp;
            for (auto& dndx : dndx_map) {
                tmp = std::max(dndx.second->Calculate(energy), 0.);
                if (dndx.first)
                    tmp /= medium.GetSumNucleons()
                        / (dndx.first->GetAtomInMolecule()
                              * dndx.first->GetAtomicNum());

                dndx_all += tmp;
            }

            return dndx_all;
        }
        else
        {
            double dndx = std::max(dndx_map[comp_ptr]->Calculate(energy), 0.);
            if (comp_ptr)
                dndx /= medium.GetSumNucleons()
                    / (comp_ptr->GetAtomInMolecule()
                          * comp_ptr->GetAtomicNum());
            return dndx;
        }
    }
    inline double CalculateStochasticLoss(
        std::shared_ptr<const Component> const& comp, double energy, double rate) override
    {
        return CalculateStochasticLoss_impl(comp, energy, rate,
            typename param_t::base_param_t::component_wise{},
            typename param_t::base_param_t::only_stochastic{});
    }
    inline size_t GetHash() const noexcept override
    {
        auto hash_digest = size_t{ 0 };
        hash_combine(hash_digest, param.GetHash(), p_def.GetHash(),
            medium.GetHash());
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
    const Medium& medium, const EnergyCutSettings& cut, Interpolant1DBuilder::Definition& def)
{
    Integral integral;
    def.function1d
        = [&integral, &param, &p_def, &medium, &cut](double energy) {
              return calculate_dedx(param, integral, p_def, medium, cut, energy,
                  typename std::decay<Param>::type::component_wise{});
          };
    def.xmin = param.GetLowerEnergyLim(p_def);
    def.rational = true;
    def.logSubst = true;
    auto hash = def.GetHash();
    hash_combine(hash, param.GetHash(), cut.GetHash());
    return Helper::InitializeInterpolation("dEdx", Interpolant1DBuilder(def), hash);
}

template <typename Param>
unique_ptr<Interpolant> build_de2dx(Param&& param, const ParticleDef& p_def,
    const Medium& medium, const EnergyCutSettings& cut, Interpolant1DBuilder::Definition& def)
{
    Integral integral;
    def.function1d
        = [&integral, &param, &p_def, &medium, &cut](double energy) {
              return calculate_de2dx(param, integral, p_def, medium, cut,
                  energy, typename std::decay<Param>::type::component_wise{});
          };
    def.xmin = param.GetLowerEnergyLim(p_def);
    auto hash = def.GetHash();
    hash_combine(hash, param.GetHash(), cut.GetHash());
    return Helper::InitializeInterpolation("dE2dx", Interpolant1DBuilder(def), hash);
}

template <typename Param, typename P, typename M>
double CrossSectionInterpolant<Param, P, M>::CalculateStochasticLoss_impl(
    std::shared_ptr<const Component> const& comp, double energy, double rate, std::true_type, std::false_type)
{
    auto weight_for_rate_in_medium = medium.GetSumNucleons()
        / (comp->GetAtomInMolecule() * comp->GetAtomicNum());
    return dndx_map[comp]->GetUpperLimit(
        energy, rate * weight_for_rate_in_medium);
}

template <typename Param, typename P, typename M>
double CrossSectionInterpolant<Param, P, M>::CalculateStochasticLoss_impl(
    std::shared_ptr<const Component> const&, double energy, double rate, std::false_type, std::false_type)
{
    return dndx_map[nullptr]->GetUpperLimit(energy, rate);
}

template <typename Param, typename P, typename M>
double CrossSectionInterpolant<Param, P, M>::CalculateStochasticLoss_impl(
    std::shared_ptr<const Component> const&, double, double, bool, std::true_type)
{
    return 1;
}

} // namespace PROPOSAL
