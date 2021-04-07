#pragma once

#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"
#include "PROPOSAL/propagation_utility/Time.h"
#include <memory>

namespace PROPOSAL {
class ExactTimeBuilder : public Time {
    std::shared_ptr<Displacement> disp;
    size_t hash;
    std::shared_ptr<UtilityIntegral> time_integral;

public:
    ExactTimeBuilder(
        std::shared_ptr<Displacement>, double mass, std::false_type);
    ExactTimeBuilder(
        std::shared_ptr<Displacement>, double mass, std::true_type);

    double FunctionToIntegral(double energy);
    double TimeElapsed(double initial_energy, double final_energy,
        double grammage, double local_density) override;
    auto GetHash() const noexcept { return hash; }
};

std::unique_ptr<Time> make_time(
    std::shared_ptr<Displacement>, const ParticleDef&, bool interpolate = true);

std::unique_ptr<Time> make_time(std::vector<std::shared_ptr<CrossSectionBase>>,
    const ParticleDef&, bool interpolate = true);

struct ApproximateTimeBuilder : public Time {
    ApproximateTimeBuilder() = default;

    template <typename... Args> ApproximateTimeBuilder(Args...) { }

    ~ApproximateTimeBuilder() = default;

    double TimeElapsed(
        double, double, double grammage, double local_density) override;
};
} // namespace PROPOSAL
