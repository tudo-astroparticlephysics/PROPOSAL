#pragma once

#include "PROPOSAL/secondaries/parametrization/photopairproduction/PhotoPairProductionInterpolant.h"
#include "PROPOSAL/crosssection/parametrization/PhotoPairProduction.h"
#include "PROPOSAL/math/Integral.h"

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
} // namespace secondaries
} // namespace PROPOSAL
