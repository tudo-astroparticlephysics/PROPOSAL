#include "PROPOSAL/math/Spherical3D.h"
#include "PROPOSAL/math/Cartesian3D.h"

using namespace PROPOSAL;

Spherical3D::Spherical3D(const nlohmann::json& config) {
    if (not config.is_array() or not(config.size() == 3))
        throw std::invalid_argument("Json array for Spherical3D is not a 3 "
                                    "component array.");

    if (!config[0].is_number() or !config[1].is_number() or !config[2].is_number())
        throw std::invalid_argument("Json array for Spherical3D must contain "
                                    "three numbers (radius, azimuth, theta).");

    for (size_t i = 0; i<3; i++)
        config[i].get_to(coordinates[i]);
}

void Spherical3D::print(std::ostream& os) const {
    os << "radius: " << coordinates[SphericalCoordinate::Radius] << "\t";
    os << "azimuth: " << coordinates[SphericalCoordinate::Azimuth] << "\t";
    os << "zenith: " << coordinates[SphericalCoordinate::Zenith] << "\n";
}

double Spherical3D::magnitude() const {
    return coordinates[SphericalCoordinate::Radius];
}

void Spherical3D::normalize() {
    coordinates[SphericalCoordinate::Radius] = 1.;
}

std::array<double, 3> Spherical3D::GetSphericalCoordinates() const {
    return coordinates;
}

std::array<double, 3> Spherical3D::GetCartesianCoordinates() const {
    auto cos_a = cos(coordinates[SphericalCoordinate::Azimuth]);
    auto sin_a = sin(coordinates[SphericalCoordinate::Azimuth]);
    auto cos_z = cos(coordinates[SphericalCoordinate::Zenith]);
    auto sin_z = sin(coordinates[SphericalCoordinate::Zenith]);
    auto r = coordinates[SphericalCoordinate::Radius];

    auto x = r * cos_a * sin_z;
    auto y = r * sin_a * sin_z;
    auto z = r * cos_z;

    return {x, y, z};
}