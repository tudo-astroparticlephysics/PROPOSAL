#pragma once

#include "PROPOSAL/secondaries/weakinteraction/WeakInteraction.h"

using std::array;
using std::vector;

namespace PROPOSAL {
namespace secondaries {
    struct NaivWeakInteraction : public secondaries::WeakInteraction {
        static constexpr int n_rnd = 0;
        const int weak_partner_type;

        NaivWeakInteraction(const ParticleDef&);

        vector<Loss::secondary_t> CalculateSecondaries(
            Loss::secondary_t, array<double, n_rnd>);
    };
} // namespace secondaries
} // namespace PROPOSAL
