#pragma once

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"

#include <memory>

namespace PROPOSAL {
class CrossSectionDEDX {
protected:
    size_t hash;
    logger_ptr logger;

public:
    template <typename T1, typename T2>
    CrossSectionDEDX(T1 const& _param, ParticleDef const& _p, T2 const& _target,
        EnergyCutSettings const& _cut)
        : hash(0)
        , logger(Logging::Get("PROPOSAL.CrossSectionDEDX"))
    {
        logger->info("Building {} {} dEdx for target {}.", _p.name, _param.name,
            _target.GetName());

        hash_combine(hash, _param.GetHash(), _p.GetHash(), _target.GetHash(),
            _cut.GetHash());
    }

    virtual ~CrossSectionDEDX() = default;

    virtual double Calculate(double energy) = 0;

    size_t GetHash() const noexcept { return hash; }
};
} // namespace PROPOSAL
