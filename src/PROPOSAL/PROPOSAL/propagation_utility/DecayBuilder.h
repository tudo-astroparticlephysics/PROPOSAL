#pragma once
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/propagation_utility/Decay.h"
#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include <cmath>
#include <memory>

namespace PROPOSAL {
class DecayBuilder : public Decay {
    using disp_ptr = std::shared_ptr<Displacement>;

    std::shared_ptr<UtilityIntegral> decay_integral;

public:
    DecayBuilder(disp_ptr, double, double, std::true_type);
    DecayBuilder(disp_ptr, double, double, std::false_type);

    double EnergyDecay(double energy, double rnd, double density) override;
};

std::unique_ptr<Decay> make_decay(
    std::shared_ptr<Displacement>, const ParticleDef&, bool interpol = true);

std::unique_ptr<Decay> make_decay(
    std::vector<std::shared_ptr<CrossSectionBase>> const&, const ParticleDef&,
    bool interpol = true);
} // namespace PROPOSAL
