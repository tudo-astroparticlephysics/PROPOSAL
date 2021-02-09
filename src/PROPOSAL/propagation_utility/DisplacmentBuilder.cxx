#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
using namespace PROPOSAL;

template <> void DisplacementBuilder<UtilityIntegral>::build_tables() { }

template <> void DisplacementBuilder<UtilityInterpolant>::build_tables()
{
    disp_integral.BuildTables("disp_", false);
}
