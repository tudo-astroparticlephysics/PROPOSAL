#pragma once

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include <string>
#include <vector>
#include <utility>

using std::string;

namespace PROPOSAL {

struct Displacement {
    static Interpolant1DBuilder::Definition interpol_def;

    Displacement() = default;
    virtual ~Displacement() = default;

    template <typename Cross> double FunctionToIntegral(Cross&&, double);
    virtual double SolveTrackIntegral(double, double) = 0;
    virtual double UpperLimitTrackIntegral(double, double) = 0;
};


template <typename Cross>
double Displacement::FunctionToIntegral(Cross&& cross, double energy)
{
    auto result = 0.0;
    for (auto& cr : cross)
        result += cr->CalculatedEdx(energy);

    return -1.0 / result;
}

struct CrossSectionVector {
    template <typename CrossVec> static double GetLowerLim(CrossVec&& cross_vec)
    {
        using val_t = typename std::decay<CrossVec>::type::value_type;
        auto result = std::max_element(
            cross_vec.begin(), cross_vec.end(), [](val_t a, val_t b) {
                return a->GetLowerEnergyLim() > b->GetLowerEnergyLim();
            });
        return (*result)->GetLowerEnergyLim();
    }

    template <typename CrossVec> static size_t GetHash(CrossVec&& cross)
    {
        auto hash_digest = size_t{ 0 };
        for (auto& c : cross)
            hash_combine(hash_digest, c->GetHash());
        return hash_digest;
    }
};

} // namespace PROPOSAL
