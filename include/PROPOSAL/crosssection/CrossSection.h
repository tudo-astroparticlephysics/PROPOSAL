
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

using Components::Component;

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
    template <typename T, typename... Args>
    inline auto build_dndx(
        std::false_type, bool interpol, T target, Args&&... args)
    {
        return std::unordered_map<std::shared_ptr<const Component>,
            std::unique_ptr<CrossSectionDNDX>> { { nullptr },
            { make_dndx(interpol, target, std::forward<Args>(args)...) } };
    }

    template <typename T, typename... Args>
    inline auto build_dndx(
        std::true_type, bool interpol, T target, Args&&... args)
    {
        auto dndx = std::unordered_map<std::shared_ptr<const Component>,
            std::unique_ptr<CrossSectionDNDX>> {};
        for (auto& c : target.GetComponents()) {
            dndx.emplace(std::make_shared<const Components::Component>(c),
                make_dndx(interpol, c, std::forward<Args>(args)...));
        }
        return dndx;
    }

    template <typename T1, typename T2, typename T3, typename... Args>
    inline auto build_dedx(std::true_type, bool interpol, T1&& param,
        T2&& p_def, T3 const& target, Args&&... args)
    {
        auto dedx = std::vector<std::unique_ptr<CrossSectionDEDX>> {};
        for (auto& c : target.GetComponents()) {
            auto calc = make_dedx(interpol, std::forward<T1>(param),
                std::forward<T2>(p_def), c, args...);
            dedx.push_back(std::move(calc));
        }
        return dedx;
    }

    template <typename T1, typename T2, typename T3, typename... Args>
    inline auto build_dedx(std::false_type, bool interpol, T1&& param,
        T2&& p_def, T3 const& target, Args&&... args)
    {
        return std::vector<std::unique_ptr<CrossSectionDEDX>> { make_dedx(
            interpol, std::forward<T1>(param), std::forward<T2>(p_def), target,
            args...) };
    }

    inline auto reweight_dndx(Medium const& m, Component const& c)
    {
        return m.GetSumNucleons() / (c.GetAtomInMolecule() * c.GetAtomicNum());
    }
}

template <class P, class M> class CrossSection : public CrossSectionBase {
public:
    P p_def;
    M medium;
    dndx_map_t dndx;
    std::vector<std::unique_ptr<CrossSectionDEDX>> dedx;

    template <typename T1, typename T2>
    CrossSection(P _p_def, M _medium, T1&& _dndx, T2&& _dedx)
        : p_def(_p_def)
        , medium(_medium)
        , dndx(std::forward<T1>(_dndx))
        , dedx(std::forward<T2>(_dedx))
    {
    }

    virtual ~CrossSection() = default;

    template <typename T>
    auto CalculatedNdx_impl(double energy, std::true_type, T comp_ptr)
    {
        if (comp_ptr)
            return dndx[comp_ptr]->Calculate(energy)
                / detail::reweight_dndx(medium, *comp_ptr);
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
        return dndx[nullptr]->Calculate(energy);
    }

    auto CalculateStochasticLoss_impl(
        std::shared_ptr<const Component> const& comp, double energy,
        double rate, std::true_type, std::false_type)
    {
        return dndx[comp]->GetUpperLimit(
            energy, rate * detail::reweight_dndx(medium, *comp));
    }
    auto CalculateStochasticLoss_impl(
        std::shared_ptr<const Component> const& comp, double energy,
        double rate, std::false_type, std::false_type)
    {
        return dndx[comp]->GetUpperLimit(energy, rate);
    }
    auto CalculateStochasticLoss_impl(std::shared_ptr<const Component> const&,
        double, double, bool, std::true_type)
    {
        return 1;
    }
    /* public: */
    /* virtual double CalculatedEdx(double) = 0; */
    /* virtual double CalculatedE2dx(double) = 0; */
    /* virtual rates_t CalculatedNdx(double) = 0; */
    /* virtual double CalculateStochasticLoss(std::shared_ptr<const Component>
     * const&, double, double) = 0; */
    /* virtual double GetLowerEnergyLim() const = 0; */

    /* virtual size_t GetHash() const noexcept = 0; */
    /* virtual InteractionType GetInteractionType() const noexcept = 0; */
};

template <typename P, typename M>
using crosssection_t
    = CrossSection<typename std::decay<P>::type, typename std::decay<M>::type>;

template <typename P, typename M>
using cross_t_ptr = std::unique_ptr<crosssection_t<P, M>>;

template <typename P, typename M>
using crosssection_list_t = std::vector<std::shared_ptr<crosssection_t<P, M>>>;

template <typename T>
size_t crosssection_hasher(size_t hash_diggest, T const& obj)
{
    hash_combine(hash_diggest, obj.GetHash());
    return hash_diggest;
}

template <typename T, typename... Args>
size_t crosssection_hasher(size_t hash_diggest, T const& obj, Args... args)
{
    hash_combine(hash_diggest, obj.GetHash());
    return crosssection_hasher(hash_diggest, args...);
}

/* template<> */
/* size_t crosssection_hasher(size_t hash_diggest, std::shared_ptr<const
 * EnergyCutSettings> const& ptr); */

} // namespace PROPOSAL
