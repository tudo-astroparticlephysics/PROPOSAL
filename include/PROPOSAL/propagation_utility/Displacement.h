#pragma once

#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/CrossSectionVector.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include <string>
#include <utility>
#include <vector>

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
} // namespace PROPOSAL
