
#include "PROPOSAL/crossection/CrossSectionInterpolant.h"
#include <cmath>

namespace PROPOSAL {
double transform_relativ_loss(double v_cut, double v_max, double v)
{
    if (v < v_cut)
        return v_max;
    if (v == v_cut)
        return v_cut;
    return v_cut * std::exp(v * std::log(v_max / v_cut));
}
}
