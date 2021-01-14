
#pragma once

#include <memory>

namespace PROPOSAL {
class CrossSectionDE2DX {
protected:
    size_t hash;

public:
    CrossSectionDE2DX() = default;
    virtual ~CrossSectionDE2DX() = default;

    virtual double Calculate(double energy) = 0;

    virtual size_t GetHash() {return hash;}
};
} // namespace PROPOSAL
