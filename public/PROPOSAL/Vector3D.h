#pragma once

#ifndef VECTOR3D_H
#define VECTOR3D_H

// #include <sstream>


namespace PROPOSAL{

class Vector3D
{
    public:
        // constructors
        Vector3D();
        Vector3D(const double x, const double y, const double z);
        Vector3D(const Vector3D &vector_3d);
        ~Vector3D();

        // operator functions and swap 
        Vector3D& operator=(const Vector3D &vector_3d);
        bool operator==(const Vector3D &vector_3d) const;
        bool operator!=(const Vector3D &vector_3d) const;
        void swap(Vector3D &vector_3d);
        friend std::ostream& operator<<(std::ostream &os, Vector3D const &vector_3d);

        //basic arithmetic
        friend Vector3D operator+ (const Vector3D &vec1, const Vector3D &vec2);
        friend Vector3D operator- (const Vector3D &vec1, const Vector3D &vec2);
        friend Vector3D operator* (const double factor1, const Vector3D &vec1);
        friend Vector3D operator* (const Vector3D &vec1, const double factor1);
        friend double scalar_product (const Vector3D &vec1, const Vector3D &vec2);
        friend Vector3D vector_product (const Vector3D &vec1, const Vector3D &vec2);
        Vector3D operator- () const;
        double magnitude () const;
        void normalise();

        // conversions to spherical and cylindrical coordinate
        void CalculateCartesianFromSpherical();
        void CalculateSphericalCoordinates();
        // void CalculateCartesianFromCylindrical();
        // void CalculateCylindricalCoordinates();

        // setter
        void SetCartesianCoordinates(const double x, const double y, const double z);
        void SetSphericalCoordinates(const double radius, const double azimuth, const double zenith);
        // void SetCylindricalCoordinates(const double radius, const double azimuth, const double height);

        // getter
        double GetX() const;
        double GetY() const;
        double GetZ() const;
        double GetRadius() const;
        double GetPhi() const;
        double GetTheta() const;

    private:
        double x_, y_, z_;
        double spheric_radius_, spheric_azimuth_, spheric_zenith_;
        // double cylindric_radius_, cylindric_azimuth_, cylindric_height_;

        double CalculateAzimuthFromCartesian();
};

}

#endif
