
#pragma once

#include "PROPOSAL/secondaries/parametrization/photoeffect/Photoeffect.h"

namespace PROPOSAL {
    namespace secondaries {
        struct PhotoeffectNoDeflection : public secondaries::Photoeffect,
                                 public DefaultSecondaries<PhotoeffectNoDeflection> {

            static constexpr int n_rnd = 0;

            PhotoeffectNoDeflection(const ParticleDef&, const Medium&) {}

            size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }
            std::vector<ParticleState> CalculateSecondaries(
                    StochasticLoss, const Component&, std::vector<double>&) final;
        };
    } // namespace secondaries
} // namespace PROPOSAL
