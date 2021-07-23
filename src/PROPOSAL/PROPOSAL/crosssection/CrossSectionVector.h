#pragma once

#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/particle/Particle.h"

#include <algorithm>
#include <type_traits>
#include <vector>

namespace PROPOSAL {
struct CrossSectionVector {
    CrossSectionVector() = delete;

    template <typename CrossVec>
    static std::vector<InteractionType> GetInteractionTypes(
        CrossVec const& cross)
    {
        auto v = std::vector<InteractionType>();
        for (auto const& c : cross)
            v.push_back(c->GetInteractionType());
        return v;
    }

    template <typename CrossVec>
    static double GetLowerLim(CrossVec const& cross)
    {
        using val_t = typename std::decay<CrossVec>::type::value_type;
        auto r = std::max_element(
            cross.begin(), cross.end(), [](val_t a_ptr, val_t b_ptr) {
                return a_ptr->GetLowerEnergyLim() > b_ptr->GetLowerEnergyLim();
            });
        return (*r)->GetLowerEnergyLim();
    }

    template <typename CrossVec>
    static double GetMinStochasticEnergy(CrossVec const& cross)
    {
        using val_t = typename std::decay<CrossVec>::type::value_type;
        auto r = std::min_element(
            cross.begin(), cross.end(), [](val_t a_ptr, val_t b_ptr) {
                return a_ptr->GetMinStochasticEnergy()
                    < b_ptr->GetMinStochasticEnergy();
            });
        return (*r)->GetMinStochasticEnergy();
    }

    template <typename CrossVec> static size_t GetHash(CrossVec const& cross)
    {
        auto hash_digest = size_t { 0 };
        for (auto const& c : cross)
            hash_combine(hash_digest, c->GetHash());
        return hash_digest;
    }
};
}
