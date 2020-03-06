
#include <cmath>
#include <iostream>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

//----------------------------------------------------------------------//
//------------------------- Constructors -------------------------------//
//----------------------------------------------------------------------//

// default constructor
Vector3D::Vector3D()
    : cartesian_(0,0,0)
    , spherical_(0,0,0)
// , cylindric_radius_  (0)
// , cylindric_azimuth_ (0)
// , cylindric_height_  (0)
{
}

Vector3D::Vector3D(const double x, const double y, const double z)
    : cartesian_(x, y, z)
    , spherical_(0,0,0)
// , cylindric_radius_  (0)
// , cylindric_azimuth_ (0)
// , cylindric_height_  (0)
{
}

// copy constructor
Vector3D::Vector3D(const Vector3D& vector_3d)
    : cartesian_(vector_3d.cartesian_)
    , spherical_(vector_3d.spherical_)
// , cylindric_radius_  (vector_3d.cylindric_radius_)
// , cylindric_azimuth_ (vector_3d.cylindric_azimuth_)
// , cylindric_height_  (vector_3d.cylindric_height_)
{
}

Vector3D::Vector3D(Vector3D&& other)
    : cartesian_(std::move(other.cartesian_))
    , spherical_(std::move(other.spherical_))
{
}

Vector3D::Vector3D(const nlohmann::json& config)
    : spherical_(0,0,0)
{
    if(not config.is_array()) throw std::invalid_argument("Vector3D is not a 3 component array.");
    if(not (config.size() == 3)) throw std::invalid_argument("Vector3D is not 3.");
    if(not config[0].is_number()) throw std::invalid_argument("x is not a number");
    if(not config[1].is_number()) throw std::invalid_argument("y is not a number");
    if(not config[2].is_number()) throw std::invalid_argument("z is not a number");

    config[0].get_to(cartesian_.x_);
    config[1].get_to(cartesian_.y_);
    config[2].get_to(cartesian_.z_);

    cartesian_.x_ *= 100; // cm
    cartesian_.y_ *= 100; // cm
    cartesian_.z_ *= 100; // cm
}

// destructor
Vector3D::~Vector3D() {}

//----------------------------------------------------------------------//
//-----------------------operator functions and swap--------------------//
//----------------------------------------------------------------------//

Vector3D& Vector3D::operator=(const Vector3D& vector_3d)
{
    if (this != &vector_3d)
    {
        Vector3D tmp(vector_3d);
        swap(tmp);
    }
    return *this;
}

bool Vector3D::operator==(const Vector3D& vector_3d) const
{
    if (cartesian_.x_ != vector_3d.cartesian_.x_)
        return false;
    else if (cartesian_.y_ != vector_3d.cartesian_.y_)
        return false;
    else if (cartesian_.z_ != vector_3d.cartesian_.z_)
        return false;
    else if (spherical_.radius_ != vector_3d.spherical_.radius_)
        return false;
    else if (spherical_.azimuth_ != vector_3d.spherical_.azimuth_)
        return false;
    else if (spherical_.zenith_ != vector_3d.spherical_.zenith_)
        return false;
    // else if (cylindric_radius_  != vector_3d.cylindric_radius_)  return false;
    // else if (cylindric_azimuth_ != vector_3d.cylindric_azimuth_) return false;
    // else if (cylindric_height_  != vector_3d.cylindric_height_)  return false;

    return true;
}

bool Vector3D::operator!=(const Vector3D& vector_3d) const
{
    return !(*this == vector_3d);
}

void Vector3D::swap(Vector3D& vector_3d)
{
    using std::swap;

    swap(cartesian_.x_, vector_3d.cartesian_.x_);
    swap(cartesian_.y_, vector_3d.cartesian_.y_);
    swap(cartesian_.z_, vector_3d.cartesian_.z_);
    swap(spherical_.radius_, vector_3d.spherical_.radius_);
    swap(spherical_.azimuth_, vector_3d.spherical_.azimuth_);
    swap(spherical_.zenith_, vector_3d.spherical_.zenith_);
    // swap(cylindric_radius_,  vector_3d.cylindric_radius_);
    // swap(cylindric_azimuth_, vector_3d.cylindric_azimuth_);
    // swap(cylindric_height_,  vector_3d.cylindric_height_);
}

namespace PROPOSAL {
std::ostream& operator<<(std::ostream& os, Vector3D const& vector_3d)
{
    std::stringstream ss;
    ss << " Vector3D (" << &vector_3d << ") ";
    os << Helper::Centered(60, ss.str()) << '\n';

    os << "Cartesian Coordinates (x[cm],y[cm],z[cm]):\n"
       << vector_3d.cartesian_.x_ << "\t" << vector_3d.cartesian_.y_ << "\t" << vector_3d.cartesian_.z_ << std::endl;
    os << "Spherical Coordinates (radius[cm],azimut[rad],zenith[rad]):\n"
       << vector_3d.spherical_.radius_ << "\t" << vector_3d.spherical_.azimuth_ << "\t" << vector_3d.spherical_.zenith_
       << std::endl;
    // os<<"\tCylindrical Coordinates
    // (radius,azimut,height):\t"<<vector_3d.cylindric_radius_<<"\t"<<vector_3d.cylindric_azimuth_<<"\t"<<vector_3d.cylindric_height_<<std::endl;

    os << Helper::Centered(60, "");
    return os;
}
} // namespace PROPOSAL

//----------------------------------------------------------------------//
//-----------------------operator basic arithmetic ---------------------//
//----------------------------------------------------------------------//

namespace PROPOSAL {

Vector3D operator+(const Vector3D& vec1, const Vector3D& vec2)
{
    Vector3D vector_sum;
    vector_sum.cartesian_.x_ = vec1.cartesian_.x_ + vec2.cartesian_.x_;
    vector_sum.cartesian_.y_ = vec1.cartesian_.y_ + vec2.cartesian_.y_;
    vector_sum.cartesian_.z_ = vec1.cartesian_.z_ + vec2.cartesian_.z_;
    return vector_sum;
}

Vector3D operator-(const Vector3D& vec1, const Vector3D& vec2)
{
    Vector3D vector_diff;
    vector_diff.cartesian_.x_ = vec1.cartesian_.x_ - vec2.cartesian_.x_;
    vector_diff.cartesian_.y_ = vec1.cartesian_.y_ - vec2.cartesian_.y_;
    vector_diff.cartesian_.z_ = vec1.cartesian_.z_ - vec2.cartesian_.z_;
    return vector_diff;
}

Vector3D operator*(const double factor1, const Vector3D& vec1)
{
    Vector3D product;
    product.cartesian_.x_ = factor1 * vec1.cartesian_.x_;
    product.cartesian_.y_ = factor1 * vec1.cartesian_.y_;
    product.cartesian_.z_ = factor1 * vec1.cartesian_.z_;
    return product;
}

Vector3D operator*(const Vector3D& vec1, const double factor1)
{
    Vector3D product;
    product.cartesian_.x_ = factor1 * vec1.cartesian_.x_;
    product.cartesian_.y_ = factor1 * vec1.cartesian_.y_;
    product.cartesian_.z_ = factor1 * vec1.cartesian_.z_;
    return product;
}

double operator*(const Vector3D& vec1, const Vector3D& vec2)
{
    return scalar_product(vec1, vec2);
}

double scalar_product(const Vector3D& vec1, const Vector3D& vec2)
{
    return vec1.cartesian_.x_ * vec2.cartesian_.x_ + vec1.cartesian_.y_ * vec2.cartesian_.y_ + vec1.cartesian_.z_ * vec2.cartesian_.z_;
}

Vector3D vector_product(const Vector3D& vec1, const Vector3D& vec2)
{
    Vector3D product;
    product.cartesian_.x_ = vec1.cartesian_.y_ * vec2.cartesian_.z_ - vec1.cartesian_.z_ * vec2.cartesian_.y_;
    product.cartesian_.y_ = vec1.cartesian_.z_ * vec2.cartesian_.x_ - vec1.cartesian_.x_ * vec2.cartesian_.z_;
    product.cartesian_.z_ = vec1.cartesian_.x_ * vec2.cartesian_.y_ - vec1.cartesian_.y_ * vec2.cartesian_.x_;
    return product;
}

} // namespace PROPOSAL

Vector3D Vector3D::operator-() const
{
    Vector3D vector_3d;
    vector_3d.cartesian_.x_ = -cartesian_.x_;
    vector_3d.cartesian_.y_ = -cartesian_.y_;
    vector_3d.cartesian_.z_ = -cartesian_.z_;
    return vector_3d;
}

double Vector3D::magnitude() const
{
    return std::sqrt(cartesian_.x_ * cartesian_.x_ + cartesian_.y_ * cartesian_.y_ + cartesian_.z_ * cartesian_.z_);
}

void Vector3D::normalise()
{
    double length   = std::sqrt(cartesian_.x_ * cartesian_.x_ + cartesian_.y_ * cartesian_.y_ + cartesian_.z_ * cartesian_.z_);
    cartesian_.x_              = cartesian_.x_ / length;
    cartesian_.y_              = cartesian_.y_ / length;
    cartesian_.z_              = cartesian_.z_ / length;
    spherical_.radius_ = 1;
}

void Vector3D::deflect(const double cosphi_deflect, const double theta_deflect)
{
    if(cosphi_deflect != 1 || theta_deflect != 0)
    {
        CalculateSphericalCoordinates();

        double sinphi_deflect = std::sqrt( std::max(0., (1. - cosphi_deflect) * (1. + cosphi_deflect) ));
        double tx = sinphi_deflect * std::cos(theta_deflect);
        double ty = sinphi_deflect * std::sin(theta_deflect);
        double tz = std::sqrt(std::max(1. - tx * tx - ty * ty, 0.));
        if(cosphi_deflect < 0. ){
            // Backward deflection
            tz = -tz;
        }

        double sinth, costh, sinph, cosph;
        sinth = std::sin(spherical_.zenith_);
        costh = std::cos(spherical_.zenith_);
        sinph = std::sin(spherical_.azimuth_);
        cosph = std::cos(spherical_.azimuth_);

        const Vector3D rotate_vector_x = Vector3D(costh * cosph, costh * sinph, -sinth);
        const Vector3D rotate_vector_y = Vector3D(-sinph, cosph, 0.);

        // Rotation towards all tree axes
        Vector3D new_direction( tz * *this + tx * rotate_vector_x + ty * rotate_vector_y );

        *this = new_direction;
    }
}

//----------------------------------------------------------------------//
//---------------Spherical and cylindrical coordinates------------------//
//----------------------------------------------------------------------//

void Vector3D::CalculateCartesianFromSpherical()
{
    cartesian_.x_ = spherical_.radius_ * std::cos(spherical_.azimuth_) * std::sin(spherical_.zenith_);
    cartesian_.y_ = spherical_.radius_ * std::sin(spherical_.azimuth_) * std::sin(spherical_.zenith_);
    cartesian_.z_ = spherical_.radius_ * std::cos(spherical_.zenith_);
}

void Vector3D::CalculateSphericalCoordinates()
{
    spherical_.radius_  = std::sqrt(cartesian_.x_ * cartesian_.x_ + cartesian_.y_ * cartesian_.y_ + cartesian_.z_ * cartesian_.z_);
    spherical_.azimuth_ = std::atan2(cartesian_.y_, cartesian_.x_);
    if (spherical_.radius_ > 0.)
    {
        spherical_.zenith_ = std::acos(cartesian_.z_ / spherical_.radius_);
    } else if (spherical_.radius_ == 0.)
    {
        // log_warn("If the radius is zero, the zenith is not defined! Zero is returned!");
        spherical_.zenith_ = 0.;
    } else
    {
        // log_fatal("The radius is negativ, which is not possible!");
    }
}

// void Vector3D::CalculateCartesianFromCylindrical()
// {
//     x_ = cylindric_radius_*std::cos(cylindric_azimuth_);
//     y_ = cylindric_radius_*std::sin(cylindric_azimuth_);
//     z_ = cylindric_height_;
// }

// void Vector3D::CalculateZylindricalCoordinates()
// {
//     cylindric_radius_  = std::sqrt(x_*x_ + y_*y_);
//     cylindric_azimuth_ = std::atan2(y_, x_);
//     cylindric_height_  = z_;
// }

Vector3D::CartesianCoordinates::CartesianCoordinates(double x, double y, double z):
    x_(x), y_(y), z_(z)
{}

Vector3D::CartesianCoordinates::CartesianCoordinates(const CartesianCoordinates& coordinates):
    x_(coordinates.x_), y_(coordinates.y_), z_(coordinates.z_)
{}

Vector3D::SphericalCoordinates::SphericalCoordinates(double radius, double azimuth, double zenith):
    radius_(radius), azimuth_(azimuth), zenith_(zenith)
{}

Vector3D::SphericalCoordinates::SphericalCoordinates(const SphericalCoordinates& coordinates):
    radius_(coordinates.radius_), azimuth_(coordinates.azimuth_), zenith_(coordinates.zenith_)
{}

//----------------------------------------------------------------------//
//---------------------private member function--------------------------//
//----------------------------------------------------------------------//
