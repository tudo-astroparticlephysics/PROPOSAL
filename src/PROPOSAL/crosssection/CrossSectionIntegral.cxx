#include "PROPOSAL/crosssection/CrossSectionIntegral.h"
#include "PROPOSAL/crosssection/parametrization/Ionization.h"
#include "PROPOSAL/math/Integral.h"

using std::get;

namespace PROPOSAL {

template <>
double calculate_dedx<crosssection::Ionization&>(crosssection::Ionization& param, Integral& integral,
    const ParticleDef& p_def, const Medium& medium,
    const EnergyCutSettings& cut, double energy, std::false_type)
{
    if (!param.dEdx_integration())
        return param.FunctionToDEdxIntegral(p_def, medium, energy, 0.);
    auto physical_lim = param.GetKinematicLimits(p_def, medium, energy);
    auto v_cut = cut.GetCut(physical_lim, energy);
    auto dEdx = [&param, &p_def, &medium, energy](double v) {
        return energy * param.FunctionToDEdxIntegral(p_def, medium, energy, v);
    };
    return integral.Integrate(
        get<crosssection::Parametrization::V_MIN>(physical_lim), v_cut, dEdx, 4);
}
} // namespace PROPOSAL
