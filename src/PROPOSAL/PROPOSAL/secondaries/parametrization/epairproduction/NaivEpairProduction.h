
#pragma once

#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/secondaries/parametrization/epairproduction/EpairProduction.h"


namespace PROPOSAL {
namespace secondaries {
    struct NaivEpairProduction
        : public secondaries::EpairProduction {
        // unordered_map<Component*, unique_ptr<Interpolant>> dndx;

        static constexpr int n_rnd = 2;
        size_t RequiredRandomNumbers() const noexcept { return n_rnd; }

        NaivEpairProduction(ParticleDef, Medium){};

        std::vector<ParticleState> CalculateSecondaries(StochasticLoss, const Component&,
                                                   std::vector<double>&) {
            auto sec = std::vector<ParticleState>{};
            return sec;
        };
    };
} // namespace secondaries
} // namespace PROPOSAL
