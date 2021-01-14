
#include "PROPOSAL/crosssection/CrossSectionDE2DX/CrossSectionDE2DXInterpolant.h"

using namespace PROPOSAL;
double CrossSectionDE2DXInterpolant::Calculate(double energy)
{
    return interpolant.evaluate(energy);
}
