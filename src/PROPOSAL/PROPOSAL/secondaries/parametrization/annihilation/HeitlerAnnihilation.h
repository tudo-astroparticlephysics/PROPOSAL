#pragma once

#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/parametrization/Annihilation.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/secondaries/parametrization/annihilation/Annihilation.h"

namespace PROPOSAL {
    namespace secondaries {
        class HeitlerAnnihilation : public secondaries::Annihilation,
                                    public DefaultSecondaries<HeitlerAnnihilation>  {

            double mass;
            static double f_term(double a1, double a2, double v);
            double f(double v, double a1, double a2, double v_min, double v_max, double rnd);

        public:
            static constexpr int n_rnd = 2;

            HeitlerAnnihilation(const ParticleDef& p, const Medium& medium) : mass(p.mass) {}

            double CalculateRho(double, double, const Component&) final;
            std::tuple<Cartesian3D, Cartesian3D> CalculateDirections(
                    const Vector3D&, double, double, double) final;
            std::tuple<double, double> CalculateEnergy(double, double) final;

            size_t RequiredRandomNumbers() const noexcept override { return n_rnd; }
            std::vector<ParticleState> CalculateSecondaries(
                    StochasticLoss, const Component&, std::vector<double>&) override;
        };
    } // namespace secondaries
} // namespace PROPOSAL
