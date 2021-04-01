#pragma once

#include "PROPOSAL/secondaries/parametrization/photopairproduction/PhotoTsai.h"

namespace PROPOSAL {
namespace secondaries {
    class PhotoTsaiForwardPeaked : public PhotoTsai,
        public DefaultSecondaries<PhotoTsaiForwardPeaked> {

    public:
        PhotoTsaiForwardPeaked() = default;
        PhotoTsaiForwardPeaked(ParticleDef p, Medium m)
            : PhotoTsai(p, m)
        {
        }

        std::tuple<Cartesian3D, Cartesian3D> CalculateDirections(
            const Vector3D&, double, double, const Component&,
            std::vector<double>) final;
    };
} // namespace secondaries
} // namespace PROPOSAL
