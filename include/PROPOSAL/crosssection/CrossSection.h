
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
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDX.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Components.h"
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
        for (auto& c : target.GetComponents())
            dndx.emplace(std::make_shared<const Components::Component>(c),
                make_dndx(interpol, c, std::forward<Args>(args)...));
        return dndx;
    }
}

template <class P, class M> struct CrossSection : public CrossSectionBase {

    CrossSection() = default;
    virtual ~CrossSection() = default;

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
