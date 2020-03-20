
#include "PROPOSAL/propagation_utility/Displacement.h"

using namespace PROPOSAL;

Displacement::Displacement(CrossSectionList cross): cross(cross) {}


double Displacement::FunctionToIntegral(double energy) {
    double result = 0.0;
    for (const auto& cr: cross)
        result += cr->CalculatedEdx(energy);

    return -1.0 / result;
}
