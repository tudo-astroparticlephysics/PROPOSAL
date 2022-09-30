#include "PROPOSAL/propagation_utility/DisplacementApproximation.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;

DisplacementApproximation::DisplacementApproximation(
        crossbase_list_t const& _cross)
        : Displacement(_cross)
{
}

double DisplacementApproximation::SolveTrackIntegral(
        double lower_lim, double upper_lim)
{
    return (upper_lim - lower_lim) * FunctionToIntegral(lower_lim);
}

double DisplacementApproximation::UpperLimitTrackIntegral(
        double lower_lim, double sum)
{
    return lower_lim + sum / FunctionToIntegral(lower_lim);
}

namespace PROPOSAL {
    std::unique_ptr<Displacement> make_displacement_approximation(
            std::vector<std::shared_ptr<CrossSectionBase>> const& cross)
    {
        return std::make_unique<DisplacementApproximation>(cross);
    }
} // namespace PROPOSAL
