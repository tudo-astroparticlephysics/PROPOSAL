#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"

using namespace PROPOSAL;

DisplacementBuilder::DisplacementBuilder(
    crossbase_list_t const& _cross, std::false_type)
    : Displacement(_cross)
    , disp_integral(std::make_unique<UtilityIntegral>(
          [this](double E) { return FunctionToIntegral(E); },
          this->GetLowerLim(), this->GetHash()))
{
}

DisplacementBuilder::DisplacementBuilder(
    crossbase_list_t const& _cross, std::true_type)
    : Displacement(_cross)
    , disp_integral(std::make_unique<UtilityInterpolant>(
          [this](double E) { return FunctionToIntegral(E); },
          this->GetLowerLim(), this->GetHash()))
{
    disp_integral->BuildTables("disp_", 500, false);
}

double DisplacementBuilder::SolveTrackIntegral(
    double lower_lim, double upper_lim)
{
    return disp_integral->Calculate(lower_lim, upper_lim);
}

double DisplacementBuilder::UpperLimitTrackIntegral(
    double lower_lim, double sum)
{
    return disp_integral->GetUpperLimit(lower_lim, sum);
}

namespace PROPOSAL {
std::unique_ptr<Displacement> make_displacement(
    std::vector<std::shared_ptr<CrossSectionBase>> const& cross,
    bool interpolate)
{
    if (interpolate)
        return std::make_unique<DisplacementBuilder>(cross, std::true_type {});
    return std::make_unique<DisplacementBuilder>(cross, std::false_type {});
}
} // namespace PROPOSAL
