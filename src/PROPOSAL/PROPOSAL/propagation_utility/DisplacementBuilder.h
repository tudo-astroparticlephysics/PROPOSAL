#pragma once
#include "PROPOSAL/propagation_utility/Displacement.h"

namespace PROPOSAL {
class UtilityIntegral;
} // namespace PROPOSAL

namespace PROPOSAL {
class DisplacementBuilder : public Displacement {
    std::unique_ptr<UtilityIntegral> disp_integral;

public:
    DisplacementBuilder(crossbase_list_t const&, std::false_type);
    DisplacementBuilder(crossbase_list_t const&, std::true_type);

    double SolveTrackIntegral(double lower_lim, double upper_lim) final;

    double UpperLimitTrackIntegral(double lower_lim, double sum) final;
};

std::unique_ptr<Displacement> make_displacement(
    std::vector<std::shared_ptr<CrossSectionBase>> const&,
    bool interpolate = false);
} // namespace PROPOSAL
