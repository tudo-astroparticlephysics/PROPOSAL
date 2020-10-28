#pragma once

#include "PROPOSAL/secondaries/bremsstrahlung/Bremsstrahlung.h"
#include "PROPOSAL/secondaries/RegisteredInDefault.h"

using std::array;
using std::vector;

namespace PROPOSAL {
namespace secondaries {
    class NaivBremsstrahlung : public secondaries::Bremsstrahlung,
                               public RegisteredInDefault<NaivBremsstrahlung> {
        static constexpr int n_rnd = 0;
        const int primary_lepton_type;

    public:
        NaivBremsstrahlung() = delete;
        NaivBremsstrahlung(ParticleDef p, Medium) :primary_lepton_type(p.particle_type) {};

        size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }
        vector<ParticleState> CalculateSecondaries(StochasticLoss, const Component&,
                                                   vector<double>&);
    };
} // namespace secondaries
} // namespace PROPOSAL
