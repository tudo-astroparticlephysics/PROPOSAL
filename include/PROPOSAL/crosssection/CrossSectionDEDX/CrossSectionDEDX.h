
#pragma once

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"

#include <memory>

namespace PROPOSAL {
struct CrossSectionDEDX {
    CrossSectionDEDX() = default;
    virtual ~CrossSectionDEDX() = default;

    virtual double Calculate(double energy) = 0;
};
} // namespace PROPOSAL
