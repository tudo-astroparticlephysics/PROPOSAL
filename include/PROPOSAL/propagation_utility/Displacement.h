#pragma once

#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/CrossSectionVector.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include <utility>
#include <vector>

namespace PROPOSAL {

class Displacement {
protected:
    using crossbase_list_t = std::vector<std::shared_ptr<CrossSectionBase>>;
    crossbase_list_t cross_list;

public:
    template <typename Cross>
    Displacement(Cross&& cross) : cross_list(std::begin(cross), std::end(cross)) {}
    virtual ~Displacement() = default;

    static Interpolant1DBuilder::Definition interpol_def;

    double FunctionToIntegral(double);
    virtual double SolveTrackIntegral(double, double) = 0;
    virtual double UpperLimitTrackIntegral(double, double) = 0;
};

} // namespace PROPOSAL
