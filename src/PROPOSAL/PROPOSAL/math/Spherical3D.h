#pragma once
#include "PROPOSAL/math/Vector3D.h"
#include <nlohmann/json.hpp>

namespace PROPOSAL  {
    class Cartesian3D;
    class Spherical3D : public Vector3D {
    public:
        Spherical3D() : Vector3D() {};
        Spherical3D(std::array<double, 3> val) : Vector3D(val) {};
        Spherical3D(double radius, double azimuth, double zenith)
            : Vector3D({radius, azimuth, zenith}) {};
        Spherical3D(const Vector3D& vec) : Spherical3D(vec.GetSphericalCoordinates()) {};
        Spherical3D(const nlohmann::json&);


        auto GetRadius() const {return coordinates[Radius];}
        auto GetAzimuth() const {return coordinates[Azimuth];}
        auto GetZenith() const {return coordinates[Zenith];}
        void SetRadius(double radius) {coordinates[Radius] = radius;}
        void SetAzimuth(double azimuth) {coordinates[Azimuth] = azimuth;}
        void SetZenith(double zenith) {coordinates[Zenith] = zenith;}

        double magnitude() const override;
        void normalize() override;
        std::array<double, 3> GetCartesianCoordinates() const override;
        std::array<double, 3> GetSphericalCoordinates() const override;

        enum SphericalCoordinate : int {
            Radius = 0,
            Azimuth = 1,
            Zenith = 2,
        };
    protected:
        void print(std::ostream&) const override;
     };
} // namespace PROPOSAL
