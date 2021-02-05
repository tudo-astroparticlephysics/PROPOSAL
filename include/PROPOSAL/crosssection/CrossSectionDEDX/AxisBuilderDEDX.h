#include "PROPOSAL/Constants.h"

#include <functional>

#include "CubicInterpolation/Axis.h"

namespace PROPOSAL {
class AxisBuilderDEDX {
    double low, up;
    size_t n;

public:
    AxisBuilderDEDX(double _low, double _up = UPPER_ENERGY_LIM_DEFAULT,
        size_t _nodes = NODES_DEDX_DEFAULT);

    void refine_definition_range(std::function<double(double)> func);

    std::unique_ptr<cubic_splines::ExpAxis<double>> Create() const;
};
} // namespace PROPOSAL
