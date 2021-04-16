#include "PROPOSAL/math/Spherical3D.h"
#include "PROPOSAL/math/Cartesian3D.h"

using namespace PROPOSAL;

Spherical3D::Spherical3D(const nlohmann::json& config) {
    if (!config.is_array() || !(config.size() == 3))
        throw std::invalid_argument("Json array for Spherical3D is not a 3 "
                                    "component array.");

    if (!config[0].is_number() || !config[1].is_number() || !config[2].is_number())
        throw std::invalid_argument("Json array for Spherical3D must contain "
                                    "three numbers (radius, azimuth, theta).");

    for (size_t i = 0; i<3; i++)
        config[i].get_to(coordinates[i]);
}

void Spherical3D::print(std::ostream& os) const {
    os << "radius: " << GetRadius() << "\t";
    os << "azimuth: " << GetAzimuth() << "\t";
    os << "zenith: " << GetZenith() << "\n";
}

double Spherical3D::magnitude() const {
    return GetRadius();
}

void Spherical3D::normalize() {
    SetRadius(1.);
}

std::array<double, 3> Spherical3D::GetSphericalCoordinates() const {
    return coordinates;
}

std::array<double, 3> Spherical3D::GetCartesianCoordinates() const {
    auto cos_a = cos(GetAzimuth());
    auto sin_a = sin(GetAzimuth());
    auto cos_z = cos(GetZenith());
    auto sin_z = sin(GetZenith());
    auto r = GetRadius();

    auto x = r * cos_a * sin_z;
    auto y = r * sin_a * sin_z;
    auto z = r * cos_z;

    return {x, y, z};
}
