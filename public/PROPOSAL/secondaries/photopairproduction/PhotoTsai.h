#pragma once

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/secondaries/DefaultFactory.h"
#include "PROPOSAL/secondaries/photopairproduction/PhotopairProduction.h"

#include "PROPOSAL/crossection/CrossSectionDNDX/CrossSectionDNDXBuilder.h"
#include "PROPOSAL/crossection/parametrization/PhotoPairProduction.h"

using PROPOSAL::Components::Component;

namespace PROPOSAL {
namespace secondaries {
    class PhotoTsai : public PhotopairProduction,
                      public RegisteredInDefault<PhotoTsai> {
        Integral integral;
        dndx_map_t dndx;

        double FunctionToIntegral(
            double energy, double x, double theta, const Component&);

    public:
        static constexpr int n_rnd = 3;

        PhotoTsai() = default;
        PhotoTsai(ParticleDef p, Medium m)
            : dndx(build_cross_section_dndx(crosssection::PhotoPairTsai(), p, m,
                  std::make_shared<EnergyCutSettings>(0.f, 1.f, false), false))
        {
        }

        double CalculateRho(double, double, const Component&) override;
        tuple<Vector3D, Vector3D> CalculateDirections(Vector3D, double, double,
            const Component&, vector<double>) override;
        tuple<double, double> CalculateEnergy(double, double, double) override;

        size_t RequiredRandomNumbers() const noexcept final { return n_rnd; }
        vector<Loss::secondary_t> CalculateSecondaries(
            double, Loss::secondary_t, const Component&, vector<double>) final;
    };
} // namespace secondaries
} // namespace PROPOSAL
