
#pragma once

#include "PROPOSAL/secondaries/parametrization/photoeffect/Photoeffect.h"

namespace PROPOSAL {
    namespace secondaries {
        struct NaivPhotoeffect : public secondaries::Photoeffect,
                                 public DefaultSecondaries<NaivPhotoeffect> {

            static constexpr int n_rnd = 0;

            NaivPhotoeffect(const ParticleDef&, const Medium&) {}

            size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }
            std::vector<ParticleState> CalculateSecondaries(
                    StochasticLoss, const Component&, std::vector<double>&) final;
        };
    } // namespace secondaries
} // namespace PROPOSAL
