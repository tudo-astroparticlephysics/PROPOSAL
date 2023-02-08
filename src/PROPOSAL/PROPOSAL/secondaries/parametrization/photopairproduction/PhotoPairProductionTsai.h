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

    class PhotoPairProductionTsaiSauter : public PhotoPairProductionTsaiForwardPeaked {
    public:
        PhotoPairProductionTsaiSauter() = default;
        PhotoPairProductionTsaiSauter(ParticleDef p, Medium m)
            : PhotoPairProductionTsaiForwardPeaked(p, m)
        {}

        std::tuple<Cartesian3D, Cartesian3D> CalculateDirections(
                const Vector3D& dir, double energy, double rho,
                const Component& comp, double rnd1, double rnd2, double rnd3) override {

            // Sampling according to EGS 5 manual, formula (2.157)

            auto energies = CalculateEnergy(energy, rho);
            double E_electron = std::get<0>(energies);
            double E_positron = std::get<1>(energies);
            double p_electron = std::sqrt(E_electron * E_electron - ME * ME);
            double p_positron = std::sqrt(E_positron * E_positron - ME * ME);

            // TODO: should this be two times the same random numbers, or two different ones?
            auto cosphi0 = (E_electron * (2 * rnd2 - 1) + p_electron) / (p_electron * (2 * rnd2 - 1) + E_electron);
            auto cosphi1 = (E_positron * (2 * rnd3 - 1) + p_positron) / (p_positron * (2 * rnd3 - 1) + E_positron);

            auto theta0 = rnd1 * 2. * PI;
            auto theta1 = std::fmod(theta0 + PI, 2. * PI);
            if (cosphi0 == -1.)
                cosphi0 *= (-1);
            if (cosphi1 == -1.)
                cosphi1 *= (-1);
            auto dir_0 = Cartesian3D(dir);
            dir_0.deflect(cosphi0, theta0);
            auto dir_1 = Cartesian3D(dir);
            dir_1.deflect(cosphi1, theta1);
            return std::make_tuple(dir_0, dir_1);
        };
    };

} // namespace secondaries
} // namespace PROPOSAL
