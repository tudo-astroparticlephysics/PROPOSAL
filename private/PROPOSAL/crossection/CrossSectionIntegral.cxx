#include "PROPOSAL/crossection/CrossSectionIntegral.h"
#include "PROPOSAL/crossection/parametrization/Ionization.h"
#include "PROPOSAL/math/Integral.h"

using std::get;

namespace PROPOSAL {

template <>
double calculate_dedx<crosssection::Ionization&>(crosssection::Ionization& param, Integral& integral,
    const ParticleDef& p_def, const Medium& medium,
    const EnergyCutSettings& cut, double energy, std::false_type, std::false_type)
{
    auto is_bhabha = param.name == "IonizBergerSeltzerBhabha";
    auto is_moller = param.name == "IonizBergerSeltzerMoller";
    if (is_bhabha || is_moller)
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
