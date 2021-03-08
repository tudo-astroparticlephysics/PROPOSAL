#pragma once

#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/secondaries/parametrization/photopairproduction/PhotopairProduction.h"

#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXBuilder.h"
#include "PROPOSAL/crosssection/parametrization/PhotoPairProduction.h"

using PROPOSAL::Components::Component;

namespace PROPOSAL {
namespace secondaries {
    class PhotoTsai : public PhotopairProduction,
                      public DefaultSecondaries<PhotoTsai> {
        Integral integral;
        Medium medium;
        dndx_map_t dndx;

        double FunctionToIntegral(
            double energy, double x, double theta, const Component&);

    public:
        static constexpr int n_rnd = 5;

        PhotoTsai() = default;
        PhotoTsai(ParticleDef p, Medium m)
            : medium(m),
              dndx(build_cross_section_dndx(crosssection::PhotoPairTsai(), p, medium,
                  std::make_shared<EnergyCutSettings>(0.f, 1.f, false), false))
        {
        }

        double CalculateRho(double, double, const Component&) override;
        std::tuple<Cartesian3D, Cartesian3D> CalculateDirections(
                const Vector3D&, double, double, const Component&,
                std::vector<double>) override;
        std::tuple<double, double> CalculateEnergy(double, double, double) override;

        size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }
        std::vector<ParticleState> CalculateSecondaries(StochasticLoss, const Component&,
                                                   std::vector<double>&) final;
    };
} // namespace secondaries
} // namespace PROPOSAL
