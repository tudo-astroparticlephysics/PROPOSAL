#pragma once

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/Cartesian3D.h"
#include "PROPOSAL/medium/Components.h"

#include <cmath>

namespace PROPOSAL {
    namespace secondaries {

        struct SauterSampling {
            std::tuple<Cartesian3D, Cartesian3D> CalculateDirections(
                    const Vector3D& dir, const Component& comp, double rnd1, double rnd2, double rnd3,
                    double E_electron, double E_positron) {

                // Sampling according to EGS 5 manual, formula (2.157)

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
    }
}