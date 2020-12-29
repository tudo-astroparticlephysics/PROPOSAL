#pragma once

#include "PROPOSAL/secondaries/parametrization/weakinteraction/WeakInteraction.h"

namespace PROPOSAL {
namespace secondaries {
    class NaivWeakInteraction
        : public secondaries::WeakInteraction,
          public DefaultSecondaries<NaivWeakInteraction> {
        const int weak_partner_type;

    public:
        static constexpr int n_rnd = 0;

        NaivWeakInteraction(const ParticleDef& p, const Medium&)
            : weak_partner_type(p.weak_partner) {}

        size_t RequiredRandomNumbers() const noexcept override { return n_rnd; }
        std::vector<ParticleState> CalculateSecondaries(StochasticLoss, const Component&,
                                                   std::vector<double>&) override;
    };
} // namespace secondaries
} // namespace PROPOSAL
