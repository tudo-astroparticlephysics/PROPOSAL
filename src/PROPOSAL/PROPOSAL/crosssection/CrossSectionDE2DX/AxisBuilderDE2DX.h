#pragma once

#include "CubicInterpolation/Axis.h"
#include "PROPOSAL/Constants.h"
#include <functional>

namespace PROPOSAL {
    class AxisBuilderDE2DX {
        double low, up;
        size_t n;

    public:
        AxisBuilderDE2DX(double _low, double _up = InterpolationSettings::UPPER_ENERGY_LIM,
                        size_t _nodes = InterpolationSettings::NODES_DE2DX);

        void refine_definition_range(
                std::function<double(double)> func, unsigned int i = 0);

        std::unique_ptr<cubic_splines::ExpAxis<double>> Create() const;
    };
} // namespace PROPOSAL
