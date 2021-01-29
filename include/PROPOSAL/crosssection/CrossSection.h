
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

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/crosssection/CrossSectionDE2DX/CrossSectionDE2DXBuilder.h"
#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDXBuilder.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXBuilder.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/Particle.h"
#include <memory>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace PROPOSAL {

using rates_t = std::unordered_map<std::shared_ptr<const Component>, double>;

struct CrossSectionBase {

    virtual ~CrossSectionBase() = default;
    virtual double CalculatedEdx(double) = 0;
    virtual double CalculatedE2dx(double) = 0;
    virtual double CalculatedNdx(
        double, std::shared_ptr<const Component> = nullptr)
        = 0;
    virtual double CalculateStochasticLoss(
        std::shared_ptr<const Component> const&, double, double)
        = 0;
    virtual double GetLowerEnergyLim() const = 0;
    virtual size_t GetHash() const noexcept = 0;
    virtual InteractionType GetInteractionType() const noexcept = 0;
    virtual std::vector<std::shared_ptr<const Component>>
    GetTargets() const noexcept = 0;
};

namespace detail {
    inline auto weight_component(Medium const& m, Component const& c)
    {
        return m.GetSumNucleons() / (c.GetAtomInMolecule() * c.GetAtomicNum());
    }

    template <typename T, typename... Args>
    inline auto build_dndx(
        std::false_type, bool interpol, T target, Args&&... args)
    {
        using comp_ptr = std::shared_ptr<const Component>;
        using dndx_ptr = std::unique_ptr<CrossSectionDNDX>;
        auto m = std::unordered_map<comp_ptr, std::tuple<double, dndx_ptr>>();
        auto calc = make_dndx(interpol, target, std::forward<Args>(args)...);
        m.emplace(comp_ptr(nullptr), std::make_tuple(1., std::move(calc)));
        return m;
    }

    template <typename T, typename... Args>
    inline auto build_dndx(
        std::true_type, bool interpol, T target, Args&&... args)
    {
        using comp_ptr = std::shared_ptr<const Component>;
        using dndx_ptr = std::unique_ptr<CrossSectionDNDX>;
        auto m = std::unordered_map<comp_ptr, std::tuple<double, dndx_ptr>>();
        for (auto& c : target.GetComponents()) {
            auto comp = std::make_shared<const Components::Component>(c);
            auto weight = weight_component(target, c);
            auto calc = make_dndx(interpol, c, std::forward<Args>(args)...);
            m.emplace(comp, std::make_tuple(weight, std::move(calc)));
        }
        return m;
    }

    template <typename Cont, typename T1, typename T2, typename T3,
        typename... Args>
    inline auto _build_dedx(Cont container, std::true_type, bool interpol,
        T1&& param, T2&& p_def, T3 const& target, Args&&... args)
    {
        for (auto& c : target.GetComponents())
            container->emplace_back(weight_component(target, c),
                make_dedx(interpol, std::forward<T1>(param),
                    std::forward<T2>(p_def), c, args...));
    }

    template <typename Cont, typename T1, typename T2, typename T3,
        typename... Args>
    inline auto _build_dedx(Cont container, std::false_type, bool interpol,
        T1&& param, T2&& p_def, T3 const& target, Args&&... args)
    {
        container->emplace_back(1.,
            make_dedx(interpol, std::forward<T1>(param),
                std::forward<T2>(p_def), target, args...));
    }

    template <typename Type, typename T1, typename T2, typename T3, typename T4,
        typename... Args>
    inline auto build_dedx(Type comp_wise, bool interpol, T1&& param,
        T2&& p_def, T3&& target, T4 cut_ptr, Args&&... args)
    {
        using dedx_ptr = std::unique_ptr<CrossSectionDEDX>;
        auto dedx = std::unique_ptr<std::vector<std::tuple<double, dedx_ptr>>>(
            nullptr);
        if (cut_ptr) {
            dedx
                = std::make_unique<std::vector<std::tuple<double, dedx_ptr>>>();
            _build_dedx(dedx.get(), comp_wise, interpol,
                std::forward<T1>(param), std::forward<T2>(p_def),
                std::forward<T3>(target), *cut_ptr,
                std::forward<Args>(args)...);
        }
        return dedx;
    }

    template <typename Cont, typename T1, typename T2, typename T3,
        typename... Args>
    inline auto _build_de2dx(Cont container, std::true_type, bool interpol,
        T1&& param, T2&& p_def, T3 const& target, Args&&... args)
    {
        for (auto& c : target.GetComponents())
            container->emplace_back(weight_component(target, c),
                make_de2dx(interpol, std::forward<T1>(param),
                    std::forward<T2>(p_def), c, args...));
    }

    template <typename Cont, typename T1, typename T2, typename T3,
        typename... Args>
    inline auto _build_de2dx(Cont container, std::false_type, bool interpol,
        T1&& param, T2&& p_def, T3 const& target, Args&&... args)
    {
        container->emplace_back(1.,
            make_de2dx(interpol, std::forward<T1>(param),
                std::forward<T2>(p_def), target, args...));
    }

    template <typename Type, typename T1, typename T2, typename T3, typename T4,
        typename... Args>
    inline auto build_de2dx(Type comp_wise, bool interpol, T1&& param,
        T2&& p_def, T3&& target, T4 cut_ptr, Args&&... args)
    {
        using de2dx_ptr = std::unique_ptr<CrossSectionDE2DX>;
        using de2dx_container = std::vector<std::tuple<double, de2dx_ptr>>;
        auto de2dx = std::unique_ptr<de2dx_container>(nullptr);
        if (!cut_ptr)
            return de2dx;
        if (cut_ptr->GetContRand()) {
            de2dx = std::make_unique<de2dx_container>();
            _build_de2dx(de2dx.get(), comp_wise, interpol,
                std::forward<T1>(param), std::forward<T2>(p_def),
                std::forward<T3>(target), *cut_ptr,
                std::forward<Args>(args)...);
        }
        return de2dx;
    }

}

template <class P, class M> class CrossSection : public CrossSectionBase {
    using comp_ptr = std::shared_ptr<const Component>;
    using dndx_ptr = std::unique_ptr<CrossSectionDNDX>;
    using dedx_ptr = std::unique_ptr<CrossSectionDEDX>;
    using de2dx_ptr = std::unique_ptr<CrossSectionDE2DX>;

protected:
    size_t hash;

    std::unordered_map<comp_ptr, std::tuple<double, dndx_ptr>> dndx;
    std::unique_ptr<std::vector<std::tuple<double, dedx_ptr>>> dedx;
    std::unique_ptr<std::vector<std::tuple<double, de2dx_ptr>>> de2dx;

    double lower_energy_lim;
    InteractionType interaction_type;

public:
    template <typename Param, typename T1, typename T2, typename T3>
    CrossSection(
        Param _param, P _p_def, M, T1&& _dndx, T2&& _dedx, T3&& _de2dx)
        : hash(0)
        , dndx(std::forward<T1>(_dndx))
        , dedx(std::forward<T2>(_dedx))
        , de2dx(std::forward<T3>(_de2dx))
        , lower_energy_lim(_param.GetLowerEnergyLim(_p_def))
        , interaction_type(static_cast<InteractionType>(
              crosssection::ParametrizationId<Param>::value))
    {
        auto logger = Logging::Get("PROPOSAL.CrossSection");
        logger->info("Building {} crosssection.",
            crosssection::ParametrizationName<Param>::value);
        for (auto const& i : dndx)
            hash_combine(hash, std::get<1>(i.second)->GetHash());
        if (dedx) {
            for (auto const& i : *dedx)
                hash_combine(hash, std::get<1>(i)->GetHash());
        }
        if (de2dx) {
            for (auto const& i : *de2dx)
                hash_combine(hash, std::get<1>(i)->GetHash());
        }
    }

    virtual ~CrossSection() = default;

protected:
    template <typename T>
    auto CalculatedNdx_impl(double energy, std::true_type, T comp_ptr)
    {
        if (comp_ptr) {
            return std::get<1>(dndx[comp_ptr])->Calculate(energy)
                / std::get<0>(dndx[comp_ptr]);
        }
        return CalculatedNdx_impl(energy, std::true_type {}, nullptr);
    }
    auto CalculatedNdx_impl(double energy, std::true_type tt, std::nullptr_t)
    {
        auto dNdx_all = 0.;
        for (auto& it : dndx)
            dNdx_all += CalculatedNdx_impl(energy, tt, it.first);
        return dNdx_all;
    }
    template <typename T>
    auto CalculatedNdx_impl(double energy, std::false_type, T)
    {
        return std::get<1>(dndx[nullptr])->Calculate(energy);
    }

    auto CalculateStochasticLoss_impl(
        std::shared_ptr<const Component> const& comp, double energy,
        double rate, std::true_type, std::false_type)
    {
        return std::get<1>(dndx[comp])
            ->GetUpperLimit(energy, rate * std::get<0>(dndx[comp]));
    }
    auto CalculateStochasticLoss_impl(
        std::shared_ptr<const Component> const& comp, double energy,
        double rate, std::false_type, std::false_type)
    {
        return std::get<1>(dndx[comp])->GetUpperLimit(energy, rate);
    }
    auto CalculateStochasticLoss_impl(std::shared_ptr<const Component> const&,
        double, double, bool, std::true_type)
    {
        return 1.;
    }

public:
    inline double CalculatedEdx(double energy) override
    {
        auto loss = 0.;
        if (not dedx)
            return loss;
        for (auto& [weight, calc] : *dedx)
            loss += calc->Calculate(energy) / weight;
        return loss;
    }

    inline double CalculatedE2dx(double energy) override
    {
        auto loss = 0.;
        if (not de2dx)
            return loss;
        for (auto& [weight, calc] : *de2dx)
            loss += calc->Calculate(energy) / weight;
        return loss;
    }

    std::vector<std::shared_ptr<const Component>>
    GetTargets() const noexcept final
    {
        std::vector<std::shared_ptr<const Component>> targets;
        for (auto& it : this->dndx)
            targets.emplace_back(it.first);
        return targets;
    }

    size_t GetHash() const noexcept final { return hash; }

    inline double GetLowerEnergyLim() const override
    {
        return lower_energy_lim;
    }

    inline InteractionType GetInteractionType() const noexcept override
    {
        return interaction_type;
    }
};

template <typename P, typename M>
using crosssection_t
    = CrossSection<typename std::decay<P>::type, typename std::decay<M>::type>;

template <typename P, typename M>
using cross_t_ptr = std::unique_ptr<crosssection_t<P, M>>;

template <typename P, typename M>
using crosssection_list_t = std::vector<std::shared_ptr<crosssection_t<P, M>>>;
} // namespace PROPOSAL
