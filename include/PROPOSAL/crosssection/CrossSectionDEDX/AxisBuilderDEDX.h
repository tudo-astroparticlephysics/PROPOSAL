#pragma once

#include "CubicInterpolation/Axis.h"
#include "PROPOSAL/Constants.h"
#include <functional>

namespace PROPOSAL {
class AxisBuilderDEDX {
    double low, up;
    size_t n;

public:
    AxisBuilderDEDX(double _low, double _up = InterpolationSettings::UPPER_ENERGY_LIM,
        size_t _nodes = InterpolationSettings::NODES_DEDX);

    void refine_definition_range(
        std::function<double(double)> func, unsigned int i = 0);

    std::unique_ptr<cubic_splines::ExpAxis<double>> Create() const;
};
} // namespace PROPOSAL
