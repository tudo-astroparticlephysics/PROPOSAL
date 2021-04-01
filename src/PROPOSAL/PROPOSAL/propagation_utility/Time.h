#pragma once

#include "PROPOSAL/math/InterpolantBuilder.h"

namespace PROPOSAL {
    class Displacement;
}

namespace PROPOSAL {
class Time {

protected:
    double mass;

public:
    Time() = default;

    Time(double);

    virtual ~Time() = default;

    virtual double TimeElapsed(double, double, double, double) = 0;
};
}
