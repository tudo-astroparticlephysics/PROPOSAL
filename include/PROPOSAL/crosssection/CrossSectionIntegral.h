
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

/* #include "PROPOSAL/EnergyCutSettings.h" */
#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/math/Integral.h"
/* #include "PROPOSAL/medium/Medium.h" */
/* #include "PROPOSAL/particle/ParticleDef.h" */

namespace PROPOSAL {
template <typename Param, typename P, typename M>
class CrossSectionIntegral : public crosssection_t<P, M> {

    using param_t = typename std::decay<Param>::type;
    using particle_t = typename std::decay<P>::type;
    using medium_t = typename std::decay<M>::type;
    using base_param_ref_t = typename std::add_lvalue_reference<
        typename param_t::base_param_t>::type;

    Integral integral;
    param_t param;
    std::shared_ptr<const EnergyCutSettings> cut;

public:
    CrossSectionIntegral(Param _param, P _p_def, M _medium,
        std::shared_ptr<const EnergyCutSettings> _cut)
        : crosssection_t<P, M>(_p_def, _medium,
            detail::build_dndx(
                typename param_t::base_param_t::component_wise {}, false,
                _medium, _param, _p_def, _cut))
        , param(_param)
        , cut(_cut)
    {
        if (typename param_t::only_stochastic {} == true and cut != nullptr) {
            throw std::invalid_argument(
                "CrossSections of parametrizations that are only stochastic do "
                "not use"
                "EnergyCuts. Pass a nullptr as an EnergyCut instead.");
        }
    }

    /* template <typename T> */
    /* auto CalculatedNdx_impl(double energy, std::true_type, T comp_ptr) */
    /* { */
    /*     if (comp_ptr) */
    /*         return dndx[comp_ptr]->Calculate(energy) */
    /*             / detail::reweight_dndx(medium, *comp_ptr); */
    /*     return CalculatedNdx_impl(energy, std::true_type {}, nullptr); */
    /* } */
    /* auto CalculatedNdx_impl(double energy, std::true_type tt, std::nullptr_t)
     */
    /* { */
    /*     auto dNdx_all = 0.; */
    /*     for (auto& it : dndx) */
    /*         dNdx_all += CalculatedNdx_impl(energy, tt, *it.first); */
    /*     return dNdx_all; */
    /* } */
    /* template <typename T> */
    /* auto CalculatedNdx_impl(double energy, std::false_type, T) */
    /* { */
    /*     return dndx[nullptr]->Calculate(energy); */
    /* } */

    /* auto CalculateStochasticLoss_impl( */
    /*     std::shared_ptr<const Component> const& comp, double energy, */
    /*     double rate, std::true_type, std::false_type) */
    /* { */
    /*     return dndx[comp]->GetUpperLimit( */
    /*         energy, rate * detail::reweight_dndx(medium, *comp)); */
    /* } */
    /* auto CalculateStochasticLoss_impl( */
    /*     std::shared_ptr<const Component> const& comp, double energy, */
    /*     double rate, std::false_type, std::false_type) */
    /* { */
    /*     return dndx[comp]->GetUpperLimit(energy, rate); */
    /* } */
    /* auto CalculateStochasticLoss_impl(std::shared_ptr<const Component>
     * const&, */
    /*     double, double, bool, std::true_type) */
    /* { */
    /*     return 1; */
    /* } */

    std::vector<std::shared_ptr<const Component>>
    GetTargets() const noexcept final
    {
        std::vector<std::shared_ptr<const Component>> targets;
        for (auto& it : this->dndx)
            targets.emplace_back(it.first);
        return targets;
    }
    inline double CalculatedEdx(double energy) override
    {
        if (cut == nullptr)
            return 0;
        auto aux = calculate_dedx(reinterpret_cast<base_param_ref_t>(param),
            integral, this->p_def, this->medium, *cut, energy,
            typename param_t::component_wise {});
        return aux;
    }
    inline double CalculatedE2dx(double energy) override
    {
        if (cut == nullptr)
            return 0;
        return calculate_de2dx(reinterpret_cast<base_param_ref_t>(param),
            integral, this->p_def, this->medium, *cut, energy,
            typename param_t::component_wise {});
    }
    inline double CalculatedNdx(double energy,
        std::shared_ptr<const Component> comp_ptr = nullptr) override
    {
        return this->CalculatedNdx_impl(energy,
            typename param_t::base_param_t::component_wise {}, comp_ptr);
    }
    inline double CalculateStochasticLoss(
        std::shared_ptr<const Component> const& comp, double energy,
        double rate) override
    {
        assert(rate > 0);
        return this->CalculateStochasticLoss_impl(comp, energy, rate,
            typename param_t::base_param_t::component_wise {},
            typename param_t::base_param_t::only_stochastic {});
    }
    inline size_t GetHash() const noexcept override
    {
        auto hash_digest = size_t { 0 };
        hash_combine(hash_digest, param.GetHash(), this->p_def.mass,
            std::abs(this->p_def.charge), this->medium.GetHash());

        // Only for WeakInteraction, the sign of the charge is important
        if (param.interaction_type == InteractionType::WeakInt)
            hash_combine(hash_digest, this->p_def.charge);

        if (cut)
            hash_combine(hash_digest, cut->GetHash());
        return hash_digest;
    }
    inline double GetLowerEnergyLim() const override
    {
        return param.GetLowerEnergyLim(this->p_def);
    }
    inline InteractionType GetInteractionType() const noexcept override
    {
        return param.interaction_type;
    }
};
} // namespace PROPOSAL

namespace PROPOSAL {
template <typename Param>
double calculate_dedx(Param&& param, Integral& integral,
    const ParticleDef& p_def, const Medium& medium,
    const EnergyCutSettings& cut, double energy, std::true_type)
{
    double sum = 0.;
    for (const auto& c : medium.GetComponents()) {
        auto lim = param.GetKinematicLimits(p_def, c, energy);
        auto v_cut = cut.GetCut(lim, energy);
        auto loss = integrate_dedx(integral, param, p_def, c, energy,
            std::get<crosssection::Parametrization::V_MIN>(lim), v_cut);
        auto weight_for_loss_in_medium = medium.GetSumNucleons()
            / (c.GetAtomInMolecule() * c.GetAtomicNum());

        sum += loss / detail::reweight_dndx(medium, c);
    }
    return energy * sum;
}

template <typename Param>
double calculate_dedx(Param&& param, Integral& integral,
    const ParticleDef& p_def, const Medium& medium,
    const EnergyCutSettings& cut, double energy, std::false_type)
{
    auto lim = param.GetKinematicLimits(p_def, medium, energy);
    auto v_cut = cut.GetCut(lim, energy);
    return integrate_dedx(integral, param, p_def, medium, energy,
               std::get<crosssection::Parametrization::V_MIN>(lim), v_cut)
        * energy;
}

namespace crosssection {
    class Ionization;
}
template <>
double calculate_dedx<crosssection::Ionization&>(crosssection::Ionization&,
    Integral&, const ParticleDef&, const Medium&, const EnergyCutSettings&,
    double, std::false_type);

template <typename Param>
double calculate_de2dx(Param&& param, Integral& integral,
    const ParticleDef& p_def, const Medium& medium,
    const EnergyCutSettings& cut, double energy, std::true_type)
{
    double sum = 0.;
    for (const auto& c : medium.GetComponents()) {
        auto lim = param.GetKinematicLimits(p_def, c, energy);
        auto v_cut = cut.GetCut(lim, energy);
        auto loss2 = integrate_de2dx(integral, param, p_def, c, energy,
            std::get<crosssection::Parametrization::V_MIN>(lim), v_cut);
        auto weight_for_loss_in_medium = medium.GetSumNucleons()
            / (c.GetAtomInMolecule() * c.GetAtomicNum());

        sum += loss2 / weight_for_loss_in_medium;
    }
    return energy * energy * sum;
}

template <typename Param>
double calculate_de2dx(Param&& param, Integral& integral,
    const ParticleDef& p_def, const Medium& medium,
    const EnergyCutSettings& cut, double energy, std::false_type)
{
    auto lim = param.GetKinematicLimits(p_def, medium, energy);
    auto v_cut = cut.GetCut(lim, energy);
    return integrate_de2dx(integral, param, p_def, medium, energy,
               std::get<crosssection::Parametrization::V_MIN>(lim), v_cut)
        * energy * energy;
}

} // namepace PROPOSAL
