#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDXInterpolant.h"
#include "PROPOSAL/Constants.h"

using namespace PROPOSAL;

std::string CrossSectionDEDXInterpolant::gen_path() const
{
    return std::string(InterpolationSettings::TABLES_PATH);
}

std::string CrossSectionDEDXInterpolant::gen_name() const
{
    return std::string("dedx_") + std::to_string(GetHash())
        + std::string(".txt");
}

double CrossSectionDEDXInterpolant::Calculate(double E) const
{
    if (E < lower_energy_lim)
        return 0.;
    return interpolant.evaluate(E);
}
