#pragma once
#include "PROPOSAL/propagation_utility/Displacement.h"

namespace PROPOSAL {
    class DisplacementApproximation : public Displacement {

    public:
        DisplacementApproximation(crossbase_list_t const&);

        double SolveTrackIntegral(double lower_lim, double upper_lim) final;

        double UpperLimitTrackIntegral(double lower_lim, double sum) final;
    };

    std::unique_ptr<Displacement> make_displacement_approximation(
            std::vector<std::shared_ptr<CrossSectionBase>> const&);
} // namespace PROPOSAL
