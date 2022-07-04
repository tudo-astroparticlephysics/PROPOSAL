#include "PROPOSAL/crosssection/CrossSectionDNDX/AxisBuilderDNDX.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/math/MathMethods.h"

using namespace PROPOSAL;

AxisBuilderDNDX::energy_limits AxisBuilderDNDX::refine_definition_range(energy_limits limits,
    std::function<double(double)> func, unsigned int i)
{
    auto ax = energy_axis_t(limits.low, limits.up, limits.nodes);
    while (!(func(ax.back_transform(i)) > 0) && i < limits.nodes)
        ++i;
    if (i == limits.nodes)
        throw std::logic_error("No positive values to build dNdx tables!");

    if (i == 0)
        limits.low = ax.back_transform(0);
    else {
        double i_accuracy = 0.1;
        auto f = [&func, &ax](double i) { return func(ax.back_transform(i)); };
        auto i_low = Bisection(f, i - 1, i, i_accuracy, 100).first;
        limits.low = ax.back_transform(i_low + i_accuracy);
    }
    return limits;
}

std::array<std::unique_ptr<AxisBuilderDNDX::axis_t>, 2>
AxisBuilderDNDX::Create(v_limits v_lim, energy_limits energy_lim)
{
    Logging::Get("CrossSection.DNDX.AxisBuilder")->trace("Create Axis with:");
    Logging::Get("CrossSection.DNDX.AxisBuilder")->trace("v_limits(low: {}, up: {} nodes: {})", v_lim.low, v_lim.up,
        v_lim.nodes);
    Logging::Get("CrossSection.DNDX.AxisBuilder")->trace("energy_limits(low: {}, up: {} nodes: {})", energy_lim.low,
        energy_lim.up, energy_lim.nodes);

    auto axis = std::array<std::unique_ptr<axis_t>, 2>();
    axis[0] = std::unique_ptr<axis_t>(std::make_unique<energy_axis_t>(
        energy_lim.low, energy_lim.up, energy_lim.nodes));
    axis[1] = std::unique_ptr<axis_t>(
        std::make_unique<v_axis_t>(v_lim.low, v_lim.up, v_lim.nodes));
    return axis;
}

std::unique_ptr<AxisBuilderDNDX::axis_t>
AxisBuilderDNDX::Create(energy_limits energy_lim)
{
    Logging::Get("CrossSection.DNDX.AxisBuilder")->trace("Create Axis with:");
    Logging::Get("CrossSection.DNDX.AxisBuilder")->trace("energy_limits(low: {}, up: {} nodes: {})", energy_lim.low,
                                                         energy_lim.up, energy_lim.nodes);

    auto axis = std::unique_ptr<axis_t>(std::make_unique<energy_axis_t>(
            energy_lim.low, energy_lim.up, energy_lim.nodes));
    return axis;
}