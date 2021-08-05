#pragma once

#include <PROPOSAL/crosssection/CrossSection.h>
#include "PROPOSAL/secondaries/parametrization/weakinteraction/WeakInteraction.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXBuilder.h"
#include "PROPOSAL/crosssection/parametrization/WeakInteraction.h"

namespace PROPOSAL {
namespace secondaries {
    class WeakCooperSarkarMertsch
        : public secondaries::WeakInteraction,
          public DefaultSecondaries<WeakCooperSarkarMertsch> {
        const int weak_partner_type;

        size_t generate_hash(const ParticleDef&, const Medium&);

        using dndx_ptr = std::unique_ptr<CrossSectionDNDX>;
        std::unique_ptr<
                std::unordered_map<size_t, std::tuple<double, dndx_ptr>>> dndx;
    public:
        static constexpr int n_rnd = 1;

        WeakCooperSarkarMertsch(const ParticleDef& p, const Medium& m)
            : weak_partner_type(p.weak_partner),
            dndx(detail::build_dndx(
                    std::true_type {}, true,
                    crosssection::WeakCooperSarkarMertsch(), p, m, nullptr,
                    generate_hash(p, m))){}

        double CalculateRelativeLoss(double, double, const Component&);
        size_t RequiredRandomNumbers() const noexcept override { return n_rnd; }
        std::vector<ParticleState> CalculateSecondaries(StochasticLoss, const Component&,
                                                   std::vector<double>&) override;
    };
} // namespace secondaries
} // namespace PROPOSAL
