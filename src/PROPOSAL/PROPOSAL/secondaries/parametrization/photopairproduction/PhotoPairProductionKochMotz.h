#pragma once

#include "PROPOSAL/secondaries/parametrization/photopairproduction/PhotoPairProductionInterpolant.h"
#include "PROPOSAL/crosssection/parametrization/PhotoPairProduction.h"

namespace PROPOSAL {
    namespace secondaries {
        class PhotoPairProductionKochMotzForwardPeaked : public PhotoPairProductionInterpolant<crosssection::PhotoPairKochMotz>,
                public DefaultSecondaries<PhotoPairProductionKochMotzForwardPeaked> {

        public:
            PhotoPairProductionKochMotzForwardPeaked() = default;
            PhotoPairProductionKochMotzForwardPeaked(ParticleDef p, Medium m)
                : PhotoPairProductionInterpolant<crosssection::PhotoPairKochMotz>(p, m)
            {}

        };
    } // namespace secondaries
} // namespace PROPOSAL
