
#include "PROPOSAL/crosssection/CrossSectionDE2DX/CrossSectionDE2DXInterpolant.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;

std::string CrossSectionDE2DXInterpolant::gen_path() const
{
    return std::string(InterpolationSettings::TABLES_PATH);
}

std::string CrossSectionDE2DXInterpolant::gen_name() const
{
    return std::string("de2dx_") + std::to_string(GetHash())
        + std::string(".txt");
}
double CrossSectionDE2DXInterpolant::Calculate(double energy) const
{
    return interpolant.evaluate(energy);
}
