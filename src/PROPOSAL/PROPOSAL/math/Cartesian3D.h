#pragma once
#include "PROPOSAL/math/Vector3D.h"
#include <nlohmann/json.hpp>

namespace PROPOSAL{
    class Spherical3D;
    class Cartesian3D : public Vector3D {
    public:
        Cartesian3D() = default;
        Cartesian3D(std::array<double, 3> val) : Vector3D(val) {};
        Cartesian3D(double x, double y, double z) : Vector3D({x, y, z}) {};
        Cartesian3D(const Vector3D& vec) : Cartesian3D(vec.GetCartesianCoordinates()) {};
        Cartesian3D(const nlohmann::json&);

        auto GetX() const {return coordinates[0];}
        auto GetY() const {return coordinates[1];}
        auto GetZ() const {return coordinates[2];}
        void SetX(double x) {coordinates[0] = x;}
        void SetY(double y) {coordinates[1] = y;}
        void SetZ(double z) {coordinates[2] = z;}

        friend Cartesian3D operator +(const Cartesian3D&, const Cartesian3D&);
        friend Cartesian3D operator -(const Cartesian3D&, const Cartesian3D&);
        friend double operator *(const Cartesian3D&, const Cartesian3D&);
        friend Cartesian3D operator*(const Cartesian3D&, double);
        friend Cartesian3D operator*(double, const Cartesian3D&);
        Cartesian3D operator-() const;
        Cartesian3D& operator+=(const Cartesian3D&);
        friend Cartesian3D vector_product(const Cartesian3D&, const Cartesian3D&);

        void deflect(double, double);
        double magnitude() const override;
        void normalize() override;
        std::array<double, 3> GetCartesianCoordinates() const override;
        std::array<double, 3> GetSphericalCoordinates() const override;

    protected:
        void print(std::ostream&) const override;
    };

} // namespace PROPOSAL
