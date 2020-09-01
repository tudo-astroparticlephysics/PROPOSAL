#pragma once

#include "PROPOSAL/secondaries/weakinteraction/WeakInteraction.h"

using std::array;
using std::vector;

namespace PROPOSAL {
namespace secondaries {
    class NaivWeakInteraction
        : public secondaries::WeakInteraction,
          public RegisteredInDefault<NaivWeakInteraction> {
        const int weak_partner_type;

    public:
        static constexpr int n_rnd = 0;

        NaivWeakInteraction(const ParticleDef& p, const Medium&)
            : weak_partner_type(p.weak_partner) {}

        size_t RequiredRandomNumbers() const noexcept { return n_rnd; }
        vector<Loss::secondary_t> CalculateSecondaries(double,
            Loss::secondary_t, const Component&, vector<double>) override;
    };
} // namespace secondaries
} // namespace PROPOSAL
