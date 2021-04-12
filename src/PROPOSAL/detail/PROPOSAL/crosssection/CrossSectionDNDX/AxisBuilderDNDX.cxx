#include "PROPOSAL/crosssection/CrossSectionDNDX/AxisBuilderDNDX.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/math/MathMethods.h"

using namespace PROPOSAL;

AxisBuilderDNDX::AxisBuilderDNDX(v_limits _v_lim, energy_limits _energy_lim)
    : logger(Logging::Get("CrossSection.DNDX.AxisBuidler"))
    , v_lim(_v_lim)
    , energy_lim(_energy_lim)
{
    logger->trace("Initialized with:");
    logger->trace("v_limits(low: {}, up: {} nodes: {})", v_lim.low, v_lim.up,
        v_lim.nodes);
    logger->trace("energy_limits(low: {}, up: {} nodes: {})", energy_lim.low,
        energy_lim.up, energy_lim.nodes);
}

void AxisBuilderDNDX::refine_definition_range(
    std::function<double(double)> func, unsigned int i)
{
    auto ax = energy_axis_t(energy_lim.low, energy_lim.up, energy_lim.nodes);
    while (!(func(ax.back_transform(i)) > 0) && i < energy_lim.nodes)
        ++i;
    if (i == energy_lim.nodes)
        throw std::logic_error("No positive values to build dNdx tables!");

    if (i == 0)
        energy_lim.low = ax.back_transform(0);
    else {
        double i_accuracy = 0.1;
        auto f = [&func, &ax](double i) { return func(ax.back_transform(i)); };
        auto i_low = Bisection(f, i - 1, i, i_accuracy, 100);
        energy_lim.low = ax.back_transform(i_low + i_accuracy);
    }
}

std::array<std::unique_ptr<AxisBuilderDNDX::axis_t>, 2>
AxisBuilderDNDX::Create() const
{
    logger->trace("Create Axis with:");
    logger->trace("v_limits(low: {}, up: {} nodes: {})", v_lim.low, v_lim.up,
        v_lim.nodes);
    logger->trace("energy_limits(low: {}, up: {} nodes: {})", energy_lim.low,
        energy_lim.up, energy_lim.nodes);

    auto axis = std::array<std::unique_ptr<axis_t>, 2>();
    axis[0] = std::unique_ptr<axis_t>(std::make_unique<energy_axis_t>(
        energy_lim.low, energy_lim.up, energy_lim.nodes));
    axis[1] = std::unique_ptr<axis_t>(
        std::make_unique<v_axis_t>(v_lim.low, v_lim.up, v_lim.nodes));
    return axis;
}
