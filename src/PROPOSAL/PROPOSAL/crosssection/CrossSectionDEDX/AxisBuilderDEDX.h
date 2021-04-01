#pragma once

#include "CubicInterpolation/Axis.h"
#include "PROPOSAL/Constants.h"

#include <exception>
#include <functional>
#include <memory>

namespace PROPOSAL {
class AxisBuilderDEDX {
    double low, up;
    size_t n;

public:
    AxisBuilderDEDX(double _low,
        double _up = InterpolationSettings::UPPER_ENERGY_LIM,
        size_t _nodes = InterpolationSettings::NODES_DEDX);

    void refine_definition_range(
        std::function<double(double)> func, unsigned int i = 0);

    std::unique_ptr<cubic_splines::ExpAxis<double>> Create() const;
};

static constexpr auto axis_builder_dedx_err_str
    = "dEdx axis builder is corrupt";
struct exception_axis_builder_dedx : public std::exception {
    virtual const char* what() const noexcept
    {
        return axis_builder_dedx_err_str;
    }
};

static constexpr auto axis_builder_dedx_out_of_range_err_str
    = "No node was found where the function value was greater than zero. "
      "Either there are no continuous losses for the cut, or something went "
      "wrong.";
struct exception_axis_builder_dedx_out_of_range
    : public exception_axis_builder_dedx {
    virtual const char* what() const noexcept
    {
        return axis_builder_dedx_out_of_range_err_str;
    }
};

} // namespace PROPOSAL
