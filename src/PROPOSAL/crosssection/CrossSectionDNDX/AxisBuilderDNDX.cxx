#include "PROPOSAL/crosssection/CrossSectionDNDX/AxisBuilderDNDX.h"

using namespace PROPOSAL;

AxisBuilderDNDX::AxisBuilderDNDX(v_limits _v_lim, energy_limits _energy_lim)
    : v_lim(_v_lim)
    , energy_lim(_energy_lim)
{
}

void AxisBuilderDNDX::refine_definition_range(
    std::function<double(double)> func, unsigned int i)
{
    auto ax = energy_axis_t(energy_lim.low, energy_lim.up, energy_lim.nodes);
    while (not(func(ax.back_transform(i)) > 0) and i < energy_lim.nodes)
        ++i;
    if (i == energy_lim.nodes)
        throw std::logic_error("No positive values to build dNdx tables!");
    energy_lim.low = ax.back_transform(i);
}

std::array<std::unique_ptr<AxisBuilderDNDX::axis_t>, 2>
AxisBuilderDNDX::Create() const
{
    auto axis = std::array<std::unique_ptr<axis_t>, 2>();
    axis[0] = std::unique_ptr<axis_t>(std::make_unique<energy_axis_t>(
        energy_lim.low, energy_lim.up, energy_lim.nodes));
    axis[1] = std::unique_ptr<axis_t>(
        std::make_unique<v_axis_t>(v_lim.low, v_lim.up, v_lim.nodes));
    return axis;
}
