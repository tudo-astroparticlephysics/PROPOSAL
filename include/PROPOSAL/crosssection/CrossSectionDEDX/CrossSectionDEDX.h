#pragma once

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"

#include <memory>

namespace PROPOSAL {
class CrossSectionDEDX {
protected:
    size_t hash;

public:
    template <typename T1, typename T2>
    CrossSectionDEDX(T1 const& param, ParticleDef const& p_def,
        T2 const& target, EnergyCutSettings const& cut)
        : hash(0)
    {
        hash_combine(hash, param.GetHash(), p_def.GetHash(), target.GetHash(),
            cut.GetHash());
    }

    virtual ~CrossSectionDEDX() = default;

    virtual double Calculate(double energy) = 0;

    size_t GetHash() const noexcept { return hash; }
};
} // namespace PROPOSAL
