#pragma once

#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/secondaries/parametrization/photopairproduction/PhotopairProduction.h"

#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXBuilder.h"
#include "PROPOSAL/crosssection/parametrization/PhotoPairProduction.h"

namespace PROPOSAL {
namespace secondaries {
    class PhotoTsai : public PhotopairProduction {

        using comp_ptr = std::shared_ptr<const Component>;
        using dndx_ptr = std::unique_ptr<CrossSectionDNDX>;

        Integral integral;
        Medium medium;
        std::unique_ptr<
            std::unordered_map<size_t, std::tuple<double, dndx_ptr>>>
            dndx;

        double FunctionToIntegral(
            double energy, double x, double theta, const Component&);

    public:
        static constexpr int n_rnd = 5;

        PhotoTsai() = default;
        PhotoTsai(ParticleDef p, Medium m)
            : medium(m)
            , dndx(detail::build_dndx(std::true_type {}, true,
                  crosssection::PhotoPairTsai(), p, medium, nullptr))
        {
        }

        double CalculateRho(double, double, const Component&) override;
        std::tuple<Cartesian3D, Cartesian3D> CalculateDirections(
            const Vector3D&, double, double, const Component&,
            std::vector<double>) override;
        std::tuple<double, double> CalculateEnergy(
            double, double, double) override;

        size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }
        std::vector<ParticleState> CalculateSecondaries(
            StochasticLoss, const Component&, std::vector<double>&) final;
    };
} // namespace secondaries
} // namespace PROPOSAL
