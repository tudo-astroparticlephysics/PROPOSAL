#pragma once

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/AxisBuilderDNDX.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXIntegral.h"
#include "PROPOSAL/methods.h"

#include <type_traits>

#include "CubicInterpolation/BicubicSplines.h"
#include "CubicInterpolation/Interpolant.h"

namespace PROPOSAL {

double transform_relative_loss(double v_cut, double v_max, double v, double c = 1);
double retransform_relative_loss(double v_cut, double v_max, double v, double c = 1);
double transform_loss_log(double v_cut, double v_max, double v);
double retransform_loss_log(double v_cut, double v_max, double v);

namespace crosssection {
    struct Compton;
    struct Ionization;
}

// general transformation rules for v nodes
// TODO: in C++17, this section can be rewritten much cleaner using if constexpr

template <typename Param, std::enable_if_t<!(std::is_base_of<crosssection::Ionization, Param>::value || std::is_base_of<crosssection::Compton, Param>::value), bool> = true>
double transform_loss(double v_cut, double v_max, double v) {
    return transform_relative_loss(v_cut, v_max, v);
}

template <typename Param, std::enable_if_t<!(std::is_base_of<crosssection::Ionization, Param>::value || std::is_base_of<crosssection::Compton, Param>::value), bool> = true>
double retransform_loss(double v_cut, double v_max, double v) {
    return retransform_relative_loss(v_cut, v_max, v);
}

// specifications for Ionization

template <typename Param, std::enable_if_t<std::is_base_of<crosssection::Ionization, Param>::value, bool> = true>
double transform_loss(double v_cut, double v_max, double v) {
    return transform_relative_loss(v_cut, v_max, v, 1.5);
}

template <typename Param, std::enable_if_t<std::is_base_of<crosssection::Ionization, Param>::value, bool> = true>
double retransform_loss(double v_cut, double v_max, double v) {
    return retransform_relative_loss(v_cut, v_max, v, 1.5);
}

// specifications for Compton

template <typename Param, std::enable_if_t<std::is_base_of<crosssection::Compton, Param>::value, bool> = true>
double transform_loss(double v_cut, double v_max, double v) {
    return transform_loss_log(v_cut, v_max, v);
}

template <typename Param, std::enable_if_t<std::is_base_of<crosssection::Compton, Param>::value, bool> = true>
double retransform_loss(double v_cut, double v_max, double v) {
    return retransform_loss_log(v_cut, v_max, v);
}

template <typename T1, typename... Args>
auto build_dndx_def(T1 const& param, ParticleDef const& p, Args... args)
{
    auto dndx = std::make_shared<CrossSectionDNDXIntegral>(param, p, args...);
    auto v_lim = AxisBuilderDNDX::v_limits { 0, 1,
        InterpolationSettings::NODES_DNDX_V };
    auto energy_lim = AxisBuilderDNDX::energy_limits();
    energy_lim.low = param.GetLowerEnergyLim(p);
    energy_lim.up = InterpolationSettings::UPPER_ENERGY_LIM;
    energy_lim.nodes = InterpolationSettings::NODES_DNDX_E ;
    auto energy_lim_refined = AxisBuilderDNDX::refine_definition_range(energy_lim,
        [dndx](double E) { return dndx->Calculate(E); });

    auto def = cubic_splines::BicubicSplines<double>::Definition();
    def.axis = AxisBuilderDNDX::Create(v_lim, energy_lim_refined);
    def.f = [dndx](double energy, double v) {
        auto lim = dndx->GetIntegrationLimits(energy);
        v = transform_loss<T1>(lim.min, lim.max, v);
        return dndx->Calculate(energy, v);
    };
    def.approx_derivates = true;
    return def;
}

class CrossSectionDNDXInterpolant : public CrossSectionDNDX, public LogTableCreation {
    std::function<double(double, double, double)> transform_v;
    std::function<double(double, double, double)> retransform_v;
    cubic_splines::Interpolant<cubic_splines::BicubicSplines<double>>
        interpolant;
    InteractionType type_id;

    std::string gen_path() const;
    std::string gen_name() const;
    size_t gen_hash(size_t) const;
    double evaluate_interpolant(double E, double vbar);

public:
    template <typename Param, typename Target>
    CrossSectionDNDXInterpolant(Param param, ParticleDef const& p,
        Target const& t, std::shared_ptr<const EnergyCutSettings> cut,
        size_t hash = 0)
        : CrossSectionDNDX(param, p, t, cut, gen_hash(hash)), LogTableCreation(gen_path(), gen_name())
        , transform_v(transform_loss<Param>)
        , retransform_v(retransform_loss<Param>)
        , interpolant(build_dndx_def(param, p, t, cut), gen_path(), gen_name())
        , type_id(static_cast<InteractionType>(
                crosssection::ParametrizationId<Param>::value))
    {
        lower_energy_lim
            = interpolant.GetDefinition().GetAxis().at(0)->GetLow();
    }

    double Calculate(double E) final;

    double Calculate(double E, double v) final;

    double GetUpperLimit(double, double) final;
};
} // namespace PROPOSAL
