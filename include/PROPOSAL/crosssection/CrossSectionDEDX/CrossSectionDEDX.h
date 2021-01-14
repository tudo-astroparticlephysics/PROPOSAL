
#pragma once

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"

#include <memory>

namespace PROPOSAL {
class CrossSectionDEDX {
protected:
    size_t hash;

public:
    CrossSectionDEDX() = default;
    virtual ~CrossSectionDEDX() = default;

    virtual double Calculate(double energy) = 0;

    virtual size_t GetHash() = 0;
};
} // namespace PROPOSAL
