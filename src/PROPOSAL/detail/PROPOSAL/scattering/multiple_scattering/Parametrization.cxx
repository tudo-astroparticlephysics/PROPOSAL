/*! \file   Scattering.cxx
 *   \brief  Source filefor the Scattering bug routines.
 *
 *   This version has a major bug and produces too small scattering angles.
 *
 *   \date   2013.08.19
 *   \author Tomasz Fuchs
 **/

#include <cmath>
#include <tuple>
#include <string>
#include <sstream>
#include <iostream>
#include <PROPOSAL/math/Spherical3D.h>


#include "PROPOSAL/methods.h"
#include "PROPOSAL/scattering/multiple_scattering/Parametrization.h"
#include "PROPOSAL/math/Cartesian3D.h"
#include "PROPOSAL/math/Spherical3D.h"

using namespace PROPOSAL::multiple_scattering;

Parametrization::Parametrization(double _mass)
    : mass(_mass)
{
}

bool Parametrization::operator==(const Parametrization& scattering) const
{
    if (mass != scattering.mass)
        return false;
    else
        return this->compare(scattering);
}

namespace PROPOSAL {
    namespace multiple_scattering {
        std::pair<Cartesian3D, Cartesian3D> ScatterInitialDirection(const Vector3D& direction, const ScatteringOffset& offsets) {
            auto& sx = offsets.sx;
            auto& sy = offsets.sy;
            auto& tx = offsets.tx;
            auto& ty = offsets.ty;
            auto sz = std::sqrt(std::max(1. - (sx * sx + sy * sy), 0.));
            auto tz = std::sqrt(std::max(1. - (tx * tx + ty * ty), 0.));

            auto direction_spherical = Spherical3D(direction);
            auto sinth = std::sin(direction_spherical.GetZenith());
            auto costh = std::cos(direction_spherical.GetZenith());
            auto sinph = std::sin(direction_spherical.GetAzimuth());
            auto cosph = std::cos(direction_spherical.GetAzimuth());

            auto rotate_vector_x = Cartesian3D(costh * cosph, costh * sinph, -sinth);
            auto rotate_vector_y = Cartesian3D(-sinph, cosph, 0.);

            // Rotation towards all tree axes
            auto direction_cartesian = Cartesian3D(direction);
            auto mean_direction = sz * direction_cartesian;
            mean_direction += sx * rotate_vector_x;
            mean_direction += sy * rotate_vector_y;

            // Rotation towards all tree axes
            auto final_direction = tz * direction_cartesian;
            final_direction += tx * rotate_vector_x;
            final_direction += ty * rotate_vector_y;

            return std::make_pair(mean_direction, final_direction);
        }
    }
}

namespace PROPOSAL {
std::ostream& operator<<(
    std::ostream& os, multiple_scattering::Parametrization const& scattering)
{
    std::stringstream ss;
    ss << " Multiple scattering (" << &scattering << ") ";
    os << Helper::Centered(60, ss.str()) << '\n';

    scattering.print(os);

    os << Helper::Centered(60, "");

    return os;
}
} // namespace PROPOSAL
