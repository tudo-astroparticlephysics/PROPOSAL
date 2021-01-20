#include "PROPOSAL/math/Cartesian3D.h"
#include "PROPOSAL/math/Spherical3D.h"
#include <cmath>

using namespace PROPOSAL;

Cartesian3D::Cartesian3D(const nlohmann::json& config) {
    if (not config.is_array() or not(config.size() == 3))
        throw std::invalid_argument("Json array for Cartesian3D is not a 3 "
                                    "component array.");

    if (!config[0].is_number() or !config[1].is_number() or !config[2].is_number())
        throw std::invalid_argument("Json array for Cartesian3D must contain "
                                    "three numbers (x, y, z).");

    for (size_t i = 0; i<3; i++)
        config[i].get_to(coordinates[i]);
}

void Cartesian3D::print(std::ostream& os) const {
    os << "x: " << coordinates[0] << "\t";
    os << "y: " << coordinates[1] << "\t";
    os << "z: " << coordinates[2] << "\n";
}

double Cartesian3D::magnitude() const {
    double sum = 0;
    for (auto& c : coordinates) {
        sum += c * c;
    }
    return std::sqrt(sum);
}

void Cartesian3D::normalize() {
    double length = magnitude();
    for (size_t i = 0; i<3; i++) {
        coordinates[i] /= length;
    }
}

namespace PROPOSAL {

    Cartesian3D operator+(const Cartesian3D &lhs, const Cartesian3D &rhs) {
        Cartesian3D sum;
        for (size_t i = 0; i < 3; i++) {
            sum[i] = lhs[i] + rhs[i];
        }
        return sum;
    }

    Cartesian3D operator-(const Cartesian3D &lhs, const Cartesian3D &rhs) {
        Cartesian3D diff;
        for (size_t i = 0; i < 3; i++) {
            diff[i] = lhs[i] - rhs[i];
        }
        return diff;
    }

    double operator*(const Cartesian3D &lhs, const Cartesian3D &rhs) {
        double sum = 0.;
        for (size_t i = 0; i < 3; i++) {
            sum += lhs[i] * rhs[i];
        }
        return sum;
    }

    Cartesian3D operator*(const Cartesian3D &lhs, double val) {
        Cartesian3D product;
        for (size_t i = 0; i < 3; i++) {
            product[i] = lhs[i] * val;
        }
        return product;
    }

    Cartesian3D operator*(double val, const Cartesian3D &rhs) {
        Cartesian3D product;
        for (size_t i = 0; i < 3; i++) {
            product[i] = rhs[i] * val;
        }
        return product;
    }

    Cartesian3D Cartesian3D::operator-() const {
        Cartesian3D negative;
        for (size_t i = 0; i < 3; i++) {
            negative[i] = -coordinates[i];
        }
        return negative;
    }

    Cartesian3D &Cartesian3D::operator+=(const Cartesian3D &rhs) {
        for (size_t i = 0; i < 3; i++) {
            coordinates[i] += rhs[i];
        }
        return *this;
    }

    Cartesian3D vector_product(const Cartesian3D& lhs, const Cartesian3D& rhs) {
        Cartesian3D product;
        product[0] = lhs[1] * rhs[2] - lhs[2] * rhs[1];
        product[1] = lhs[2] * rhs[0] - lhs[0] * rhs[2];
        product[2] = lhs[0] * rhs[1] - lhs[1] * rhs[0];
        return product;
    }
}

void Cartesian3D::deflect(double cosphi_deflect, double theta_deflect) {
    if(cosphi_deflect != 1 || theta_deflect != 0)
    {
        double sinphi_deflect = std::sqrt( std::max(0., (1. - cosphi_deflect) * (1. + cosphi_deflect) ));
        double tx = sinphi_deflect * std::cos(theta_deflect);
        double ty = sinphi_deflect * std::sin(theta_deflect);
        double tz = std::sqrt(std::max(1. - tx * tx - ty * ty, 0.));
        if(cosphi_deflect < 0. ){
            // Backward deflection
            tz = -tz;
        }

        double sinth, costh, sinph, cosph;
        auto spherical = this->GetSphericalCoordinates();
        sinth = std::sin(spherical[Spherical3D::SphericalCoordinate::Zenith]);
        costh = std::cos(spherical[Spherical3D::SphericalCoordinate::Zenith]);
        sinph = std::sin(spherical[Spherical3D::SphericalCoordinate::Azimuth]);
        cosph = std::cos(spherical[Spherical3D::SphericalCoordinate::Azimuth]);

        auto rotate_vector_x = Cartesian3D(costh * cosph, costh * sinph, -sinth);
        auto rotate_vector_y = Cartesian3D(-sinph, cosph, 0.);

        // Rotation towards all tree axes
        for (size_t i = 0; i < 3; i++) {
            coordinates[i] = tz * coordinates[i] + tx * rotate_vector_x[i] + ty * rotate_vector_y[i];
        }
    }
}

std::array<double, 3> Cartesian3D::GetCartesianCoordinates() const {
    return coordinates;
}

std::array<double, 3> Cartesian3D::GetSphericalCoordinates() const {
    auto r = magnitude();

    if (r == 0)
        return {0., 0., 0.};

    auto azimuth = std::atan2(coordinates[1], coordinates[0]);
    auto zenith = std::acos(coordinates[2] / r);

    return {r, azimuth, zenith};
}

