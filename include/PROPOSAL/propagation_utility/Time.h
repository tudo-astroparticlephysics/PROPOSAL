#pragma once

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/propagation_utility/Displacement.h"

namespace PROPOSAL {
class Time {
    std::shared_ptr<Displacement> disp;

protected:
    double mass;
    size_t hash;

public:
    Time() = default;

    Time(std::shared_ptr<Displacement> _disp, double _mass)
        : disp(_disp)
        , mass(_mass)
        , hash(0)
    {
        hash_combine(hash, mass);
    }
    virtual ~Time() = default;

    static Interpolant1DBuilder::Definition interpol_def;

    virtual double TimeElapsed(double, double, double, double) = 0;

    double FunctionToIntegral(double energy);

    auto GetHash() const noexcept { return hash; }
};
}
