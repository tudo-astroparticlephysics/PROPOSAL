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

    std::shared_ptr<spdlog::logger> logger;

public:
    struct v_limits {
        double low, up;
        size_t nodes;
    };

    struct energy_limits {
        double low, up;
        size_t nodes;
    };

private:
    v_limits v_lim;
    energy_limits energy_lim;

public:
    AxisBuilderDNDX(v_limits, energy_limits);

    void refine_definition_range(
        std::function<double(double)> func, unsigned int i = 0);

    std::array<std::unique_ptr<axis_t>, 2> Create() const;
};
