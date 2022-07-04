#pragma once

#include "CubicInterpolation/Axis.h"
#include "PROPOSAL/Constants.h"
#include <array>
#include <functional>
#include <spdlog/fwd.h>

class AxisBuilderDNDX {
    using axis_t = cubic_splines::Axis<double>;
    using energy_axis_t = cubic_splines::ExpAxis<double>;
    using v_axis_t = cubic_splines::LinAxis<double>;

public:
    struct v_limits {
        double low, up;
        size_t nodes;
    };

    struct energy_limits {
        double low, up;
        size_t nodes;
    };

public:
    static energy_limits refine_definition_range(energy_limits limits,
                                          std::function<double(double)> func,
                                          unsigned int i = 0);

    static std::array<std::unique_ptr<axis_t>, 2> Create(v_limits v_lim, energy_limits energy_lim);
    static std::unique_ptr<axis_t> Create(energy_limits energy_lim);

};
