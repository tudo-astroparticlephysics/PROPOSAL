#include "PROPOSAL/propagation_utility/TimeBuilder.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"

using namespace PROPOSAL;

ExactTimeBuilder::ExactTimeBuilder(
    std::shared_ptr<Displacement> _disp, double _mass, std::false_type)
    : Time(_mass)
    , disp(_disp)
    , hash(disp->GetHash())
    , time_integral(std::make_unique<UtilityIntegral>(
          [this](double E) { return FunctionToIntegral(E); },
          _disp->GetLowerLim(), this->GetHash()))
{
}

ExactTimeBuilder::ExactTimeBuilder(
    std::shared_ptr<Displacement> _disp, double _mass, std::true_type)
    : Time(_mass)
    , disp(_disp)
    , hash(disp->GetHash())
    , time_integral(std::make_unique<UtilityInterpolant>(
          [this](double E) { return FunctionToIntegral(E); },
          _disp->GetLowerLim(), this->GetHash()))
{
    time_integral->BuildTables("time_", 500, false);
}

double ExactTimeBuilder::TimeElapsed(double initial_energy, double final_energy,
    double grammage, double local_density)
{
    (void)grammage;
    assert(initial_energy >= final_energy);
    return time_integral->Calculate(initial_energy, final_energy)
        / local_density;
}

double ExactTimeBuilder::FunctionToIntegral(double energy)
{
    auto square_momentum = std::max((energy - mass) * (energy + mass), 0.);
    auto particle_momentum = std::sqrt(square_momentum);
    auto aux = disp->FunctionToIntegral(energy);
    if (aux == 0)
        return 0;
    return aux * energy / (particle_momentum * SPEED);
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
