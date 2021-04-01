#include "PROPOSAL/crosssection/CrossSectionDEDX/CrossSectionDEDXInterpolant.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/methods.h"

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

size_t CrossSectionDEDXInterpolant::gen_hash(size_t hash) const {
    hash_combine(hash,
                 InterpolationSettings::NODES_DEDX,
                 InterpolationSettings::UPPER_ENERGY_LIM);
    return hash;
}

double CrossSectionDEDXInterpolant::Calculate(double E) const
{
    if (E < lower_energy_lim)
        return 0.;
    return interpolant.evaluate(E);
}
