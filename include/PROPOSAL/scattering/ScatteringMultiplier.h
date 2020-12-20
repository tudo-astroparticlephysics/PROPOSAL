#pragma once

#include "PROPOSAL/scattering/Scattering.h"

namespace PROPOSAL {
class ScatteringMultiplier : public Scattering {
    double multiple_scatt = 1;
    std::vector<std::pair<InteractionType, double>> stochastic_deflect;

    std::array<double, 2> _scale_deflect(std::array<double, 2>& angles, InteractionType t) override
    {
        for (auto m : stochastic_deflect) {
            if (m.first == t) {
                for (auto& a : angles)
                    a *= m.second;
                return angles;
            }
        }
        return angles;
    }

    std::array<double, 4> _scale_scatter(std::array<double, 4>& angles) override
    {
        for (auto& a : angles)
            a *= multiple_scatt;
        return angles;
    }

public:
    template <typename T1, typename T2>
    ScatteringMultiplier(T1&& _m, T2&& _s, double _mm,
        std::vector<std::pair<InteractionType, double>> _sm)
        : Scattering(std::forward<T1>(_m), std::forward<T2>(_s))
        , multiple_scatt(_mm)
        , stochastic_deflect(_sm)
    {
    }
};
} // namespace PROPOSAL
