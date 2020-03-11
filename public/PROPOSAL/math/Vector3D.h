
/******************************************************************************
 *                                                                            *
 * This file is part of the simulation tool PROPOSAL.                         *
 *                                                                            *
 * Copyright (C) 2017 TU Dortmund University, Department of Physics,          *
 *                    Chair Experimental Physics 5b                           *
 *                                                                            *
 * This software may be modified and distributed under the terms of a         *
 * modified GNU Lesser General Public Licence version 3 (LGPL),               *
 * copied verbatim in the file "LICENSE".                                     *
 *                                                                            *
 * Modifcations to the LGPL License:                                          *
 *                                                                            *
 *      1. The user shall acknowledge the use of PROPOSAL by citing the       *
 *         following reference:                                               *
 *                                                                            *
 *         J.H. Koehne et al.  Comput.Phys.Commun. 184 (2013) 2070-2090 DOI:  *
 *         10.1016/j.cpc.2013.04.001                                          *
 *                                                                            *
 *      2. The user should report any bugs/errors or improvments to the       *
 *         current maintainer of PROPOSAL or open an issue on the             *
 *         GitHub webpage                                                     *
 *                                                                            *
 *         "https://github.com/tudo-astroparticlephysics/PROPOSAL"            *
 *                                                                            *
 ******************************************************************************/


#pragma once

#include <sstream>
#include "PROPOSAL/json.hpp"

namespace PROPOSAL {

class Vector3D
{
public:
    // constructors
    Vector3D();
    Vector3D(const double x, const double y, const double z);
    Vector3D(const Vector3D& vector_3d);
    Vector3D(Vector3D&& other);
    Vector3D(const nlohmann::json&);
    ~Vector3D();

    //-------------------------------------//
    // operator functions and swap
    Vector3D& operator=(const Vector3D& vector_3d);
    bool operator==(const Vector3D& vector_3d) const;
    bool operator!=(const Vector3D& vector_3d) const;
    void swap(Vector3D& vector_3d);
    friend std::ostream& operator<<(std::ostream& os, Vector3D const& vector_3d);

    //-------------------------------------//
    // basic arithmetic
    friend Vector3D operator+(const Vector3D& vec1, const Vector3D& vec2);
    friend Vector3D operator-(const Vector3D& vec1, const Vector3D& vec2);
    friend Vector3D operator*(const double factor1, const Vector3D& vec1);
    friend Vector3D operator*(const Vector3D& vec1, const double factor1);
    friend double operator*(const Vector3D& vec1, const Vector3D& vec2);
    friend double scalar_product(const Vector3D& vec1, const Vector3D& vec2);
    friend Vector3D vector_product(const Vector3D& vec1, const Vector3D& vec2);
    Vector3D operator-() const;
    double magnitude() const;
    void normalise();
    void deflect(const double , const double);

    struct CartesianCoordinates {
        CartesianCoordinates() {};
        CartesianCoordinates(double, double, double);
        CartesianCoordinates(const CartesianCoordinates&);

        double x_, y_, z_;
    };

    struct SphericalCoordinates {
        SphericalCoordinates() {};
        SphericalCoordinates(double, double, double);
        SphericalCoordinates(const SphericalCoordinates&);

        double radius_, azimuth_, zenith_;
    };

    //-------------------------------------//
    // conversions to spherical and cylindrical coordinate
    void CalculateCartesianFromSpherical();
    void CalculateSphericalCoordinates();
    // void CalculateCartesianFromCylindrical();
    // void CalculateCylindricalCoordinates();

    //-------------------------------------//
    // setter
    void SetCartesianCoordinates(const double x, const double y, const double z)
    {
        cartesian_.x_ = x;
        cartesian_.y_ = y;
        cartesian_.z_ = z;
    }
    void SetSphericalCoordinates(const double radius, const double azimuth, const double zenith)
    {
        spherical_.radius_  = radius;
        spherical_.azimuth_ = azimuth;
        spherical_.zenith_  = zenith;
    }
    // void SetCylindricalCoordinates(const double radius, const double azimuth, const double height)
    // {
    //     cylindric_radius_  = radius;
    //     cylindric_azimuth_ = azimuth;
    //     cylindric_height_  = height;
    // }

    //-------------------------------------//
    // getter
    double GetX() const { return cartesian_.x_; }
    double GetY() const { return cartesian_.y_; }
    double GetZ() const { return cartesian_.z_; }
    double GetRadius() const { return spherical_.radius_; }
    double GetPhi() const { return spherical_.azimuth_; }
    double GetTheta() const { return spherical_.zenith_; }
    Vector3D::CartesianCoordinates GetCartesianCoordinates() const { return cartesian_; }
    Vector3D::SphericalCoordinates GetSphericalCoordinates() const { return spherical_; }

    //----------------------------------------------//
private:
    CartesianCoordinates cartesian_;
    SphericalCoordinates spherical_;

    // double cylindric_radius_, cylindric_azimuth_, cylindric_height_;
};

} // namespace PROPOSAL
