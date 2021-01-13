
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

#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/CrossSectionInterpolantBase.h"
#include "PROPOSAL/math/Interpolant.h"

namespace PROPOSAL {

template <typename Param, typename P, typename M>
class CrossSectionInterpolant : public crosssection_t<P, M>,
                                public CrossSectionInterpolantBase {

    using param_t = typename std::decay<Param>::type;
    using base_param_ref_t = typename std::add_lvalue_reference<
        typename param_t::base_param_t>::type;

    param_t param;
    std::shared_ptr<const EnergyCutSettings> cut;

protected:
    std::unique_ptr<Interpolant> de2dx;

public:
    CrossSectionInterpolant(Param _param, P _p_def, M _medium,
        std::shared_ptr<const EnergyCutSettings> _cut)
        : crosssection_t<P, M>(_p_def, _medium,
            detail::build_dndx(
                typename param_t::base_param_t::component_wise {}, true, _medium,
                _param, _p_def, _cut),
            detail::build_dedx(
                typename param_t::base_param_t::component_wise {}, true,
                _param, _p_def,_medium, *_cut))
        , param(_param) // only for back transformation
        , cut(_cut)
    {
        if (typename param_t::only_stochastic {} == true and cut != nullptr) {
            throw std::invalid_argument(
                "CrossSections of parametrizations that are only "
                "stochastic do "
                "not use"
                "EnergyCuts. Pass a nullptr as an EnergyCut instead.");
        }

        if (cut != nullptr) {
            if (cut->GetContRand())
                de2dx = build_de2dx(reinterpret_cast<base_param_ref_t>(param),
                    this->p_def, this->medium, *cut, dE2dx_def);
        }
    }
    std::vector<std::shared_ptr<const Component>>
    GetTargets() const noexcept final
    {
        std::vector<std::shared_ptr<const Component>> targets;
        for (auto& it : this->dndx) {
            targets.emplace_back(it.first);
        }
        return targets;
    }
    inline double CalculatedNdx(double energy,
        std::shared_ptr<const Component> comp_ptr = nullptr) override
    {
        return this->CalculatedNdx_impl(energy,
            typename param_t::base_param_t::component_wise {}, comp_ptr);
    }

    inline double CalculatedE2dx(double energy) override
    {
        if (de2dx == nullptr)
            return 0;
        return de2dx->Interpolate(energy);
    }
    inline double CalculateStochasticLoss(
        std::shared_ptr<const Component> const& comp, double energy,
        double rate) override
    {
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

        if (cut != nullptr)
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

template <typename Param>
std::unique_ptr<Interpolant> build_de2dx(Param&& param,
    const ParticleDef& p_def, const Medium& medium,
    const EnergyCutSettings& cut, Interpolant1DBuilder::Definition& def)
{
    Integral integral;
    def.function1d = [&integral, &param, &p_def, &medium, &cut](double energy) {
        return calculate_de2dx(param, integral, p_def, medium, cut, energy,
            typename std::decay<Param>::type::component_wise {});
    };
    def.xmin = param.GetLowerEnergyLim(p_def);
    auto hash = def.GetHash();
    hash_combine(hash, param.GetHash(), p_def.mass, std::abs(p_def.charge),
        medium.GetHash(), cut.GetHash());
    if (param.interaction_type == InteractionType::WeakInt)
        hash_combine(hash, p_def.charge);
    return Helper::InitializeInterpolation(
        "dE2dx", Interpolant1DBuilder(def), hash);
}

} // namespace PROPOSAL
