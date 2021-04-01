#pragma once

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/crosssection/CrossSection.h"

#include <algorithm>
#include <type_traits>
#include <vector>

namespace PROPOSAL {
struct CrossSectionVector {
    CrossSectionVector() = delete;

    template <typename CrossVec>
    static std::vector<InteractionType> GetInteractionTypes(
        CrossVec&& cross_vec)
    {
        auto v = std::vector<InteractionType>();
        for (auto& c_ptr : cross_vec)
            v.push_back(c_ptr->GetInteractionType());
        return v;
    }

    template <typename CrossVec> static double GetLowerLim(CrossVec&& cross_vec)
    {
        using val_t = typename std::decay<CrossVec>::type::value_type;
        auto result = std::max_element(
            cross_vec.begin(), cross_vec.end(), [](val_t a_ptr, val_t b_ptr) {
                return a_ptr->GetLowerEnergyLim() > b_ptr->GetLowerEnergyLim();
            });
        return (*result)->GetLowerEnergyLim();
    }

    template <typename CrossVec> static size_t GetHash(CrossVec&& cross)
    {
        auto hash_digest = size_t{ 0 };
        for (auto& c_ptr : cross)
            hash_combine(hash_digest, c_ptr->GetHash());
        return hash_digest;
    }
};
}
