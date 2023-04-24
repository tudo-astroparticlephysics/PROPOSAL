#pragma once

#include "PROPOSAL/secondaries/parametrization/photopairproduction/PhotoPairProductionInterpolant.h"
#include "PROPOSAL/crosssection/parametrization/PhotoPairProduction.h"
#include "PROPOSAL/secondaries/parametrization/photopairproduction/SauterSampling.h"

namespace PROPOSAL {
    namespace secondaries {
        class PhotoPairProductionKochMotzForwardPeaked : public PhotoPairProductionInterpolant<crosssection::PhotoPairKochMotz> {

        public:
            PhotoPairProductionKochMotzForwardPeaked() = default;
            PhotoPairProductionKochMotzForwardPeaked(ParticleDef p, Medium m)
                : PhotoPairProductionInterpolant<crosssection::PhotoPairKochMotz>(p, m)
            {}

        };

        class PhotoPairProductionKochMotzSauter : public PhotoPairProductionKochMotzForwardPeaked, public SauterSampling,
                                                  public DefaultSecondaries<PhotoPairProductionKochMotzSauter> {
        public:
            PhotoPairProductionKochMotzSauter() = default;
            PhotoPairProductionKochMotzSauter(ParticleDef p, Medium m)
                : PhotoPairProductionKochMotzForwardPeaked(p, m)
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
