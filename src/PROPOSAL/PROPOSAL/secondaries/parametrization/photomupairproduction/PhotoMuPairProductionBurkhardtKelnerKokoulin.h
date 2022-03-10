#pragma once

#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/secondaries/parametrization/photomupairproduction/PhotoMuPairProduction.h"

#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXBuilder.h"
#include "PROPOSAL/crosssection/parametrization/PhotoMuPairProduction.h"

#include <cmath>
#include <tuple>

using std::make_tuple;

namespace PROPOSAL {
namespace secondaries {
    class PhotoMuPairProductionBurkhardtKelnerKokoulin : public PhotoMuPairProduction,
            public DefaultSecondaries<PhotoMuPairProductionBurkhardtKelnerKokoulin> {

                using dndx_ptr = std::unique_ptr<CrossSectionDNDX>;

                Medium medium;
                std::unique_ptr<std::unordered_map<size_t, std::tuple<double, dndx_ptr>>> dndx;

            public:
                static constexpr int n_rnd = 2;

                PhotoMuPairProductionBurkhardtKelnerKokoulin() = default;
                PhotoMuPairProductionBurkhardtKelnerKokoulin(ParticleDef p, Medium m)
                : medium(m)
                , dndx(detail::build_dndx(
                        std::true_type {}, true, crosssection::PhotoMuPairBurkhardtKelnerKokoulin(), p, medium, nullptr))
                                          {}

                double Calculatex(double energy, double rnd, const Component& comp) override;

                std::tuple<Cartesian3D, Cartesian3D> CalculateDirections(
                        const Vector3D& dir, double energy,
                        const Component& comp, double rnd1) override;

                std::tuple<double, double> CalculateEnergy(double energy, double x) override {
                    return make_tuple((1. - x) * energy, x * energy);
                };

                size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }

                std::vector<ParticleState> CalculateSecondaries(
                        StochasticLoss loss, const Component& comp,
                        std::vector<double>& rnd) final;
    };
} // namespace secondaries
} // namespace PROPOSAL
