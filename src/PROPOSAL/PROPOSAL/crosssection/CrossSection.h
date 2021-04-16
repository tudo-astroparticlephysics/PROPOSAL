
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
#include "PROPOSAL/crosssection/parametrization/Parametrization.hpp"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include <memory>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace PROPOSAL {

struct CrossSectionBase {

    virtual ~CrossSectionBase() = default;
    virtual double CalculatedEdx(double) = 0;
    virtual double CalculatedE2dx(double) = 0;
    virtual double CalculatedNdx(double) = 0;
    virtual double CalculatedNdx(double, size_t) = 0;
    virtual double CalculateCumulativeCrosssection(double, size_t, double) = 0;
    virtual std::vector<std::pair<size_t, double>> CalculatedNdx_PerTarget(
        double)
        = 0;
    virtual double CalculateStochasticLoss(size_t, double, double) = 0;
    virtual double GetLowerEnergyLim() const = 0;
    virtual size_t GetHash() const noexcept = 0;
    virtual InteractionType GetInteractionType() const noexcept = 0;
    virtual std::string GetParametrizationName() const noexcept = 0;
};

namespace detail {
    inline auto weight_component(Medium const& m, Component const& c)
    {
        return m.GetSumNucleons() / (c.GetAtomInMolecule() * c.GetAtomicNum());
    }

    template <typename Param>
    inline auto build_dndx(std::false_type, bool interpol, Param param,
        ParticleDef p, Medium m, std::shared_ptr<const EnergyCutSettings> cut,
        size_t hash = 0)
    {
        using dndx_ptr_t = std::unique_ptr<CrossSectionDNDX>;
        using dndx_map_t
            = std::unordered_map<size_t, std::tuple<double, dndx_ptr_t>>;
        if (cut)
            if (cut->GetEcut() == INF && cut->GetVcut() == 1)
                return std::unique_ptr<dndx_map_t>();
        auto calc = make_dndx(interpol, param, p, m, cut, hash);
        auto dndx_map = std::make_unique<dndx_map_t>();
        dndx_map->emplace(m.GetHash(), std::make_tuple(1., std::move(calc)));
        return dndx_map;
    }

    template <typename Param>
    inline auto build_dndx(std::true_type, bool interpol, Param param,
        ParticleDef p, Medium m, std::shared_ptr<const EnergyCutSettings> cut,
        size_t hash = 0)
    {
        using dndx_ptr_t = std::unique_ptr<CrossSectionDNDX>;
        using dndx_map_t
            = std::unordered_map<size_t, std::tuple<double, dndx_ptr_t>>;
        if (cut) // TODO: is this branch realy necessary, why is a dndx created
                 // for these settings?
            if (cut->GetEcut() == INF && cut->GetVcut() == 1)
                return std::unique_ptr<dndx_map_t>();
        auto dndx_map = std::make_unique<dndx_map_t>();
        for (auto& c : m.GetComponents()) {
            auto comp_hash = c.GetHash();
            auto weight = weight_component(m, c);
            auto calc = make_dndx(interpol, param, p, c, cut, hash);
            dndx_map->emplace(
                comp_hash, std::make_tuple(weight, std::move(calc)));
        }
        return dndx_map;
    }

    template <typename Cont, typename T1, typename T2, typename T3,
        typename... Args>
    inline auto _build_dedx(Cont container, std::true_type, bool interpol,
        T1&& param, T2&& p_def, T3 const& target, Args&&... args)
    {
        for (auto& c : target.GetComponents()) {
            auto weight_comp = weight_component(target, c);
            auto dedx = make_dedx(interpol, std::forward<T1>(param),
                std::forward<T2>(p_def), c, args...);
            if (dedx)
                container->emplace_back(weight_comp, std::move(dedx));
        }
    }

    template <typename Cont, typename T1, typename T2, typename T3,
        typename... Args>
    inline auto _build_dedx(Cont container, std::false_type, bool interpol,
        T1&& param, T2&& p_def, T3 const& target, Args&&... args)
    {
        auto weight_comp = 1.;
        auto dedx = make_dedx(interpol, std::forward<T1>(param),
            std::forward<T2>(p_def), target, args...);
        if (dedx)
            container->emplace_back(weight_comp, std::move(dedx));
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

    size_t generate_cross_hash(size_t, std::string, unsigned int,
        crosssection::Parametrization<Medium> const&, ParticleDef const&,
        Medium const&, std::shared_ptr<const EnergyCutSettings>);

    size_t generate_cross_hash(size_t, std::string, unsigned int,
        crosssection::Parametrization<Component> const&, ParticleDef const&,
        Medium const&, std::shared_ptr<const EnergyCutSettings>);

    double calculate_lower_energy_lim(
        std::vector<std::tuple<double, std::unique_ptr<CrossSectionDEDX>>>*);

    std::shared_ptr<spdlog::logger> init_logger(std::string const&, size_t,
        ParticleDef const&, Medium const&,
        std::shared_ptr<const EnergyCutSettings>);
}

template <typename comp_wise, typename only_stochastic>
class CrossSection : public CrossSectionBase {

    using comp_ptr = std::shared_ptr<const Component>;
    using dndx_ptr = std::unique_ptr<CrossSectionDNDX>;
    using dedx_ptr = std::unique_ptr<CrossSectionDEDX>;
    using de2dx_ptr = std::unique_ptr<CrossSectionDE2DX>;

protected:
    size_t hash;
    std::shared_ptr<spdlog::logger> logger;

    std::unique_ptr<std::unordered_map<size_t, std::tuple<double, dndx_ptr>>>
        dndx;
    std::unique_ptr<std::vector<std::tuple<double, dedx_ptr>>> dedx;
    std::unique_ptr<std::vector<std::tuple<double, de2dx_ptr>>> de2dx;

    double lower_energy_lim;
    InteractionType interaction_type;
    std::string param_name;

public:
    template <typename Param,
        typename _name = crosssection::ParametrizationName<Param>,
        typename _id = crosssection::ParametrizationId<Param>>
    CrossSection(Param param, ParticleDef p, Medium m,
        std::shared_ptr<const EnergyCutSettings> cut, bool interpol,
        size_t _hash = 0)
        : hash(detail::generate_cross_hash(
            _hash, _name::value, _id::value, param, p, m, cut))
        , logger(detail::init_logger(_name::value, _id::value, p, m, cut))
        , dndx(detail::build_dndx(
              comp_wise {}, interpol, param, p, m, cut, hash))
        , dedx(detail::build_dedx(
              comp_wise {}, interpol, param, p, m, cut, hash))
        , de2dx(detail::build_de2dx(
              comp_wise {}, interpol, param, p, m, cut, hash))
        , lower_energy_lim(param.GetLowerEnergyLim(p))
        , interaction_type(static_cast<InteractionType>(_id::value))
        , param_name(_name::value)
    {
        // initialize hash
        hash = 0;
        if (dndx) {
            for (auto& dndx_: *dndx)
                hash_combine(hash, std::get<1>(dndx_.second)->GetHash());
        }
        if (dedx) {
            for (auto& dedx_: *dedx)
                hash_combine(hash, std::get<1>(dedx_)->GetHash());
        }
        if (de2dx) {
            for (auto& de2dx_: *de2dx)
                hash_combine(hash, std::get<1>(de2dx_)->GetHash());
        }
    }

    virtual ~CrossSection() = default;

protected:
    double CalculateStochasticLoss_impl(
        size_t target_hash, double E, double rate, std::false_type)
    {
        return std::get<1>((*dndx)[target_hash])
            ->GetUpperLimit(E, rate * std::get<0>((*dndx)[target_hash]));
    }

    double CalculateStochasticLoss_impl(size_t, double, double, std::true_type)
    {
        return 1.;
    }

public:
    double CalculatedNdx(double E) override
    {
        auto dNdx_all = 0.;
        if (dndx)
            for (auto& it : *dndx) {
                dNdx_all += std::get<1>((*dndx)[it.first])->Calculate(E)
                    / std::get<0>((*dndx)[it.first]);
            }
        return dNdx_all;
    };

    double CalculatedNdx(double E, size_t target_hash) override
    {
        if (dndx)
            return std::get<1>((*dndx)[target_hash])->Calculate(E)
                / std::get<0>((*dndx)[target_hash]);
        return 0.;
    };

    double CalculateCumulativeCrosssection(
        double E, size_t hash, double v) override
    {
        if (dndx)
            return std::get<1>((*dndx)[hash])->Calculate(E, v)
                / std::get<0>((*dndx)[hash]);
        return 0.;
    }

    std::vector<std::pair<size_t, double>> CalculatedNdx_PerTarget(
        double E) override
    {
        std::vector<std::pair<size_t, double>> rates = {};
        if (dndx) {
            for (auto& c : *dndx)
                rates.push_back({ c.first, CalculatedNdx(E, c.first) });
        }
        return rates;
    }

    double CalculateStochasticLoss(size_t hash, double E, double rate) override
    {
        if (dndx)
            return CalculateStochasticLoss_impl(
                hash, E, rate, only_stochastic {});
        throw std::logic_error("Can not calculate stochastic loss if dndx "
                               "calculator is not defined. The crosssection "
                               "is probably defined to be only-continuous.");
    }

    inline double CalculatedEdx(double energy) override
    {
        auto loss = 0.;
        if (!dedx)
            return loss;
        // will produce no working in cpp17
        // for (auto& [weight, calc] : *dedx)
        //     loss += calc->Calculate(energy) / weight;
        for (auto& weight_calc : *dedx)
            loss += std::get<1>(weight_calc)->Calculate(energy)
                / std::get<0>(weight_calc);
        return loss;
    }

    inline double CalculatedE2dx(double energy) override
    {
        auto loss = 0.;
        if (!de2dx)
            return loss;
        // will produce no working in cpp17
        // for (auto& [weight, calc] : *de2dx)
        //     loss += calc->Calculate(energy) / weight;
        for (auto& weight_calc : *de2dx)
            loss += std::get<1>(weight_calc)->Calculate(energy)
                / std::get<0>(weight_calc);
        return loss;
    }

    size_t GetHash() const noexcept override { return hash; }

    inline double GetLowerEnergyLim() const override
    {
        return lower_energy_lim;
    }

    inline InteractionType GetInteractionType() const noexcept override
    {
        return interaction_type;
    }

    inline std::string GetParametrizationName() const noexcept override
    {
        return param_name;
    }
};

using cross_t_ptr = std::unique_ptr<CrossSectionBase>;

using crosssection_list_t = std::vector<std::shared_ptr<CrossSectionBase>>;
} // namespace PROPOSAL
