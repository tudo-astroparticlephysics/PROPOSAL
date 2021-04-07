#include "PROPOSAL/propagation_utility/DecayBuilder.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"

using namespace PROPOSAL;

DecayBuilder::DecayBuilder(
    disp_ptr _disp, double _lifetime, double _mass, std::false_type)
    : Decay(_disp, _lifetime, _mass)
    , decay_integral(std::make_unique<UtilityIntegral>(
          [this](double E) { return FunctionToIntegral(E); },
          disp->GetLowerLim(), this->GetHash()))
{
}

DecayBuilder::DecayBuilder(
    disp_ptr _disp, double _lifetime, double _mass, std::true_type)
    : Decay(_disp, _lifetime, _mass)
    , decay_integral(std::make_unique<UtilityInterpolant>(
          [this](double E) { return FunctionToIntegral(E); },
          disp->GetLowerLim(), this->GetHash()))
{
    decay_integral->BuildTables("decay_", 500, true);
}

double DecayBuilder::EnergyDecay(double energy, double rnd, double density)
{
    auto rndd = -std::log(rnd) * density;
    auto rnddMin
        = decay_integral->Calculate(energy, disp->GetLowerLim()) / lifetime;
    if (rndd >= rnddMin)
        return disp->GetLowerLim();
    return decay_integral->GetUpperLimit(energy, rndd * lifetime);
}

namespace PROPOSAL {
std::unique_ptr<Decay> make_decay(
    std::shared_ptr<Displacement> disp, const ParticleDef& p, bool interpol)
{
    auto lifetime = p.lifetime;
    auto mass = p.mass;
    if (interpol)
        return std::make_unique<DecayBuilder>(
            disp, lifetime, mass, std::true_type {});
    return std::make_unique<DecayBuilder>(
        disp, lifetime, mass, std::false_type {});
}

std::unique_ptr<Decay> make_decay(
    std::vector<std::shared_ptr<CrossSectionBase>> const& cross,
    const ParticleDef& p, bool interpol)
{
    auto disp = std::shared_ptr<Displacement>(make_displacement(cross, false));
    return make_decay(disp, p, interpol);
}
} // namespace PROPOSAL
