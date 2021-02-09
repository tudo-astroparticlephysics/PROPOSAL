#include "PROPOSAL/crosssection/CrossSectionDNDX/AxisBuilderDNDX.h"
#include "PROPOSAL/Logging.h"

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
    while (not(func(ax.back_transform(i+1)) > 0) and i < energy_lim.nodes)
        ++i;
    if (i == energy_lim.nodes)
        throw std::logic_error("No positive values to build dNdx tables!");
    energy_lim.low = ax.back_transform(i);
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
