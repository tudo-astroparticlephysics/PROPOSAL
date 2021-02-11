#include "PROPOSAL/propagation_utility/TimeBuilder.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"

using namespace PROPOSAL;

ExactTimeBuilder::ExactTimeBuilder(
    std::shared_ptr<Displacement> _disp, double _mass, std::false_type)
    : Time(_disp, _mass)
    , time_integral(std::make_unique<UtilityIntegral>(
          [this](double E) { return FunctionToIntegral(E); },
          _disp->GetLowerLim(), this->GetHash()))
{
}

ExactTimeBuilder::ExactTimeBuilder(
    std::shared_ptr<Displacement> _disp, double _mass, std::true_type)
    : Time(_disp, _mass)
    , time_integral(std::make_unique<UtilityInterpolant>(
          [this](double E) { return FunctionToIntegral(E); },
          _disp->GetLowerLim(), this->GetHash()))
{
    time_integral->BuildTables("time_", false);
}

double ExactTimeBuilder::TimeElapsed(double initial_energy, double final_energy,
    double grammage, double local_density)
{
    (void)grammage;
    assert(initial_energy >= final_energy);
    return time_integral->Calculate(initial_energy, final_energy)
        / local_density;
}

namespace PROPOSAL {
std::unique_ptr<Time> make_time(
    std::shared_ptr<Displacement> disp, const ParticleDef& p, bool interpol)
{
    auto time = std::unique_ptr<Time>();
    if (interpol)
        time = std::make_unique<ExactTimeBuilder>(
            disp, p.mass, std::true_type {});
    else
        time = std::make_unique<ExactTimeBuilder>(
            disp, p.mass, std::false_type {});
    return time;
}

std::unique_ptr<Time> make_time(
    std::vector<std::shared_ptr<CrossSectionBase>> cross, const ParticleDef& p,
    bool interpol)
{
    auto disp = std::shared_ptr<Displacement>(make_displacement(cross, false));
    return make_time(disp, p, interpol);
}
} // namespace PROPOSAL

double ApproximateTimeBuilder::TimeElapsed(
    double, double, double grammage, double local_density)
{
    assert(grammage >= 0);
    return grammage / (local_density * SPEED);
}
