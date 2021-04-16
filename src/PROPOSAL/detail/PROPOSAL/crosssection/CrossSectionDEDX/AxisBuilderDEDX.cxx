#include "PROPOSAL/crosssection/CrossSectionDEDX/AxisBuilderDEDX.h"
#include "PROPOSAL/math/MathMethods.h"

using namespace PROPOSAL;

AxisBuilderDEDX::AxisBuilderDEDX(double _low, double _up, size_t _node)
    : low(_low)
    , up(_up)
    , n(_node)
{
}

void AxisBuilderDEDX::refine_definition_range(
    std::function<double(double)> func, unsigned int i)
{
    auto ax = cubic_splines::ExpAxis<double>(low, up, n);
    while (!(func(ax.back_transform(i)) > 0) && i < n)
        ++i;
    if (i == n)
        throw exception_axis_builder_dedx_out_of_range();

    if (i == 0)
        return;

    double i_accuracy = 0.1;
    auto f = [&func, &ax](double i) { return func(ax.back_transform(i)); };
    auto i_low = Bisection(f, i - 1, i, i_accuracy, 100);
    low = ax.back_transform(i_low + i_accuracy);
}

std::unique_ptr<cubic_splines::ExpAxis<double>> AxisBuilderDEDX::Create() const
{
    return std::make_unique<cubic_splines::ExpAxis<double>>(low, up, n);
}
