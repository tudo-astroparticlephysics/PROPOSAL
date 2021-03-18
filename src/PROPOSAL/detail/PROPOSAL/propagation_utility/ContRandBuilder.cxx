#include "PROPOSAL/propagation_utility/ContRandBuilder.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
using namespace PROPOSAL;

namespace PROPOSAL {
template <> void ContRandBuilder<UtilityIntegral>::build_tables() { }

template <> void ContRandBuilder<UtilityInterpolant>::build_tables()
{
    cont_rand_integral.BuildTables("cont_rand_", 200, false);
}
} // namespace PROPOSAL

namespace PROPOSAL {
std::unique_ptr<ContRand> make_contrand(std::shared_ptr<Displacement> disp,
    std::vector<std::shared_ptr<CrossSectionBase>> const& cross,
    bool interpolate)
{
    auto cont_rand = std::unique_ptr<ContRand>();
    if (interpolate)
        cont_rand = std::make_unique<ContRandBuilder<UtilityInterpolant>>(
            disp, cross);
    else
        cont_rand
            = std::make_unique<ContRandBuilder<UtilityIntegral>>(disp, cross);
    return cont_rand;
}

std::unique_ptr<ContRand> make_contrand(
    std::vector<std::shared_ptr<CrossSectionBase>> const& cross,
    bool interpolate)
{
    return make_contrand(make_displacement(cross, false), cross, interpolate);
}
} // namespace PROPOSAL
