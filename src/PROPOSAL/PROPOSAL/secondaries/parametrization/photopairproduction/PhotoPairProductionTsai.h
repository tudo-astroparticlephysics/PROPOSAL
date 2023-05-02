#pragma once

#include "PROPOSAL/secondaries/parametrization/photopairproduction/PhotoPairProductionInterpolant.h"
#include "PROPOSAL/crosssection/parametrization/PhotoPairProduction.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/secondaries/parametrization/photopairproduction/SauterSampling.h"

namespace PROPOSAL {
namespace secondaries {
    class PhotoPairProductionTsaiForwardPeaked : public PhotoPairProductionInterpolant<crosssection::PhotoPairTsai> {

    public:
        PhotoPairProductionTsaiForwardPeaked() = default;
        PhotoPairProductionTsaiForwardPeaked(ParticleDef p, Medium m)
            : PhotoPairProductionInterpolant<crosssection::PhotoPairTsai>(p, m)
            {}

    };

    class PhotoPairProductionTsai : public PhotoPairProductionTsaiForwardPeaked {
    Integral integral;
    double FunctionToIntegral(double energy, double x, double theta,
                              const Component&);

    public:
        PhotoPairProductionTsai() = default;
        PhotoPairProductionTsai(ParticleDef p, Medium m)
            : PhotoPairProductionTsaiForwardPeaked(p, m) {}

        std::tuple<Cartesian3D, Cartesian3D> CalculateDirections(
                const Vector3D&, double, double, const Component&,
                double, double, double) override;
    };

    class PhotoPairProductionTsaiSauter : public PhotoPairProductionTsaiForwardPeaked, public SauterSampling {
    public:
        PhotoPairProductionTsaiSauter() = default;
        PhotoPairProductionTsaiSauter(ParticleDef p, Medium m)
            : PhotoPairProductionTsaiForwardPeaked(p, m)
        {}

        std::tuple<Cartesian3D, Cartesian3D> CalculateDirections(
                const Vector3D& dir, double energy, double rho,
                const Component& comp, double rnd1, double rnd2, double rnd3) override {
            auto energies = CalculateEnergy(energy, rho);
            double E_electron = std::get<0>(energies);
            double E_positron = std::get<1>(energies);
            return SauterSampling::CalculateDirections(dir, comp, rnd1, rnd2, rnd3, E_electron, E_positron);
        };
    };

} // namespace secondaries
} // namespace PROPOSAL
