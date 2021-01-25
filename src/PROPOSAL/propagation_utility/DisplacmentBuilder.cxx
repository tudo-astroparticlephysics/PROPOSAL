#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"

using namespace PROPOSAL;

template <> void DisplacementBuilder<UtilityIntegral>::build_tables() { }

template <> void DisplacementBuilder<UtilityInterpolant>::build_tables()
{
    auto def = cubic_splines::CubicSplines<double>::Definition();
    disp_integral.BuildTables("disp", std::move(def), false);
}
