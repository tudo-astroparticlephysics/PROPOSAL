#pragma once

#include <cstddef>

namespace PROPOSAL {
class CrossSectionDE2DX {
protected:
    size_t hash;

public:
    CrossSectionDE2DX(size_t hash);

    virtual ~CrossSectionDE2DX() = default;

    virtual double Calculate(double energy) const = 0;

    virtual size_t GetHash() const noexcept { return hash; }
};

} // namespace PROPOSAL
