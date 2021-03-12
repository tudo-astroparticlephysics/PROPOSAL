
#include "PROPOSAL/crosssection/CrossSectionDE2DX/CrossSectionDE2DXInterpolant.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/methods.h"

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

size_t CrossSectionDE2DXInterpolant::gen_hash(size_t hash) const {
    hash_combine(hash,
                 InterpolationSettings::NODES_DE2DX,
                 InterpolationSettings::UPPER_ENERGY_LIM);
    return hash;
}
double CrossSectionDE2DXInterpolant::Calculate(double energy) const
{
    if (energy < lower_energy_lim)
        return 0.;
    return interpolant.evaluate(energy);
}
