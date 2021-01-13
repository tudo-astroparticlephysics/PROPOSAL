
#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDXInterpolant.h"

using namespace PROPOSAL;
double CrossSectionDEDXInterpolant::Calculate(double energy)
{
    return interpolant.evaluate(energy);
}
