#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/propagation_utility/Displacement.h"

using namespace PROPOSAL;

double Displacement::FunctionToIntegral(double energy)
{
    auto result = 0.0;
    for (auto& cr : cross_list)
        result += cr->CalculatedEdx(energy);

    return (result > 0) ? -1.0 / result : 0.;
}
