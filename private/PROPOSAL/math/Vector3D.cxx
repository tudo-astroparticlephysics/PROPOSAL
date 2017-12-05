
#include <cmath>
#include <iostream>

#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;


//----------------------------------------------------------------------//
//------------------------- Constructors -------------------------------//
//----------------------------------------------------------------------//

// default constructor
Vector3D::Vector3D()
    : x_ (0)
    , y_ (0)
    , z_ (0)
    , spheric_radius_  (0)
    , spheric_azimuth_ (0)
    , spheric_zenith_  (0)
    // , cylindric_radius_  (0)
    // , cylindric_azimuth_ (0)
    // , cylindric_height_  (0)
{
}


Vector3D::Vector3D(const double x, const double y, const double z)
    : x_ (x)
    , y_ (y)
    , z_ (z)
    , spheric_radius_  (0)
    , spheric_azimuth_ (0)
    , spheric_zenith_  (0)
    // , cylindric_radius_  (0)
    // , cylindric_azimuth_ (0)
    // , cylindric_height_  (0)
{
}

// copy constructor
Vector3D::Vector3D(const Vector3D &vector_3d)
    : x_ (vector_3d.x_)
    , y_ (vector_3d.y_)
    , z_ (vector_3d.z_)
    , spheric_radius_  (vector_3d.spheric_radius_)
    , spheric_azimuth_ (vector_3d.spheric_azimuth_)
    , spheric_zenith_  (vector_3d.spheric_zenith_)
    // , cylindric_radius_  (vector_3d.cylindric_radius_)
    // , cylindric_azimuth_ (vector_3d.cylindric_azimuth_)
    // , cylindric_height_  (vector_3d.cylindric_height_)
{
}

// destructor
Vector3D::~Vector3D(){}


//----------------------------------------------------------------------//
//-----------------------operator functions and swap--------------------//
//----------------------------------------------------------------------//


Vector3D& Vector3D::operator=(const Vector3D &vector_3d)
{
    if (this != &vector_3d)
    {
      Vector3D tmp(vector_3d);
      swap(tmp);
    }
    return *this;
}

bool Vector3D::operator==(const Vector3D &vector_3d) const
{
    if      (x_ != vector_3d.x_) return false;
    else if (y_ != vector_3d.y_) return false;
    else if (z_ != vector_3d.z_) return false;
    else if (spheric_radius_  != vector_3d.spheric_radius_)  return false;
    else if (spheric_azimuth_ != vector_3d.spheric_azimuth_) return false;
    else if (spheric_zenith_  != vector_3d.spheric_zenith_)  return false;
    // else if (cylindric_radius_  != vector_3d.cylindric_radius_)  return false;
    // else if (cylindric_azimuth_ != vector_3d.cylindric_azimuth_) return false;
    // else if (cylindric_height_  != vector_3d.cylindric_height_)  return false;

    return true;
}

bool Vector3D::operator!=(const Vector3D &vector_3d) const
{
    return !(*this == vector_3d);
}

void Vector3D::swap(Vector3D& vector_3d)
{
    using std::swap;

    swap(x_, vector_3d.x_);
    swap(y_, vector_3d.y_);
    swap(z_, vector_3d.z_);
    swap(spheric_radius_,  vector_3d.spheric_radius_);
    swap(spheric_azimuth_, vector_3d.spheric_azimuth_);
    swap(spheric_zenith_,  vector_3d.spheric_zenith_);
    // swap(cylindric_radius_,  vector_3d.cylindric_radius_);
    // swap(cylindric_azimuth_, vector_3d.cylindric_azimuth_);
    // swap(cylindric_height_,  vector_3d.cylindric_height_);
}

namespace PROPOSAL
{
std::ostream& operator<<(std::ostream& os, Vector3D const &vector_3d)
{
    std::stringstream ss;
    ss << " Vector3D (" << &vector_3d << ") ";
    os << Helper::Centered(60, ss.str()) << '\n';

    os << "Cartesian Coordinates (x[cm],y[cm],z[cm]):\n"<<vector_3d.x_<<"\t"<<vector_3d.y_<<"\t"<<vector_3d.z_<<std::endl;
    os << "Spherical Coordinates (radius[cm],azimut[rad],zenith[rad]):\n"<<vector_3d.spheric_radius_<<"\t"<<vector_3d.spheric_azimuth_<<"\t"<<vector_3d.spheric_zenith_<<std::endl;
    // os<<"\tCylindrical Coordinates (radius,azimut,height):\t"<<vector_3d.cylindric_radius_<<"\t"<<vector_3d.cylindric_azimuth_<<"\t"<<vector_3d.cylindric_height_<<std::endl;

    os << Helper::Centered(60, "");
    return os;
}
} // namespace PROPOSAL


//----------------------------------------------------------------------//
//-----------------------operator basic arithmetic ---------------------//
//----------------------------------------------------------------------//

namespace PROPOSAL
{

Vector3D operator+(const Vector3D &vec1, const Vector3D &vec2)
{
    Vector3D vector_sum;
    vector_sum.x_ = vec1.x_ + vec2.x_;
    vector_sum.y_ = vec1.y_ + vec2.y_;
    vector_sum.z_ = vec1.z_ + vec2.z_;
    return vector_sum;
}

Vector3D operator-(const Vector3D &vec1, const Vector3D &vec2)
{
    Vector3D vector_diff;
    vector_diff.x_ = vec1.x_ - vec2.x_;
    vector_diff.y_ = vec1.y_ - vec2.y_;
    vector_diff.z_ = vec1.z_ - vec2.z_;
    return vector_diff;
}

Vector3D operator*(const double factor1, const Vector3D &vec1)
{
    Vector3D product;
    product.x_ = factor1*vec1.x_;
    product.y_ = factor1*vec1.y_;
    product.z_ = factor1*vec1.z_;
    return product;
}

Vector3D operator*(const Vector3D &vec1, const double factor1)
{
    Vector3D product;
    product.x_ = factor1*vec1.x_;
    product.y_ = factor1*vec1.y_;
    product.z_ = factor1*vec1.z_;
    return product;
}

double operator*(const Vector3D &vec1, const Vector3D &vec2)
{
    return scalar_product(vec1, vec2);
}

double scalar_product(const Vector3D &vec1, const Vector3D &vec2)
{
    return vec1.x_*vec2.x_ + vec1.y_*vec2.y_ + vec1.z_*vec2.z_;
}

Vector3D vector_product(const Vector3D &vec1, const Vector3D &vec2)
{
    Vector3D product;
    product.x_ = vec1.y_*vec2.z_ - vec1.z_*vec2.y_;
    product.y_ = vec1.z_*vec2.x_ - vec1.x_*vec2.z_;
    product.z_ = vec1.x_*vec2.y_ - vec1.y_*vec2.x_;
    return product;
}

} // namespace PROPOSAL

Vector3D Vector3D::operator- () const
{
    Vector3D vector_3d;
    vector_3d.x_ = -x_;
    vector_3d.y_ = -y_;
    vector_3d.z_ = -z_;
    return vector_3d;
}

double Vector3D::magnitude() const
{
    return std::sqrt(x_*x_ + y_*y_ + z_*z_);
}

void Vector3D::normalise()
{
    double length = std::sqrt(x_*x_ + y_*y_ + z_*z_);
    x_ = x_/length;
    y_ = y_/length;
    z_ = z_/length;
}


//----------------------------------------------------------------------//
//---------------Spherical and cylindrical coordinates------------------//
//----------------------------------------------------------------------//


void Vector3D::CalculateCartesianFromSpherical()
{
    using namespace std;
    x_ = spheric_radius_*cos(spheric_azimuth_)*sin(spheric_zenith_);
    y_ = spheric_radius_*sin(spheric_azimuth_)*sin(spheric_zenith_);
    z_ = spheric_radius_*cos(spheric_zenith_);
}

void Vector3D::CalculateSphericalCoordinates()
{
    spheric_radius_  = std::sqrt(x_*x_ + y_*y_ + z_*z_);
    spheric_azimuth_ = CalculateAzimuthFromCartesian();
    // spheric_azimuth_ = std::atan2(y_, x_);
    if (spheric_radius_ > 0.)
    {
        spheric_zenith_ = std::acos(z_/spheric_radius_);
    }
    else if (spheric_radius_ == 0.)
    {
        // log_warn("If the radius is zero, the zenith is not defined! Zero is returned!");
        spheric_zenith_ = 0.;
    }
    else
    {
        // log_fatal("The radius is negativ, which is not possible!");
    }
}

// void Vector3D::CalculateCartesianFromCylindrical()
// {
//     using namespace std;
//     x_ = cylindric_radius_*cos(cylindric_azimuth_);
//     y_ = cylindric_radius_*sin(cylindric_azimuth_);
//     z_ = cylindric_height_;
// }

// void Vector3D::CalculateZylindricalCoordinates()
// {
//     cylindric_radius_  = std::sqrt(x_*x_ + y_*y_);
//     cylindric_azimuth_ = CalculateAzimuthFromCartesian();
// //     cylindric_azimuth_ = std::atan2(y_, x_);
//     cylindric_height_  = z_;
// }


//----------------------------------------------------------------------//
//---------------------private member function--------------------------//
//----------------------------------------------------------------------//

// in principal its the std::atan2(y_, x_) function
double Vector3D::CalculateAzimuthFromCartesian()
{
    using namespace std;
    if (x_ > 0.)
    {
        return atan(y_/x_);
    }
    else if (x_ < 0.)
    {
        if (y_ >= 0.)
        {
            return atan(y_/x_) + PI;
        }
        else
        {
            return atan(y_/x_) - PI;
        }
    }
    else if (x_ == 0)
    {
        if (y_ > 0)
        {
            return PI/2.;
        }
        else if (y_< 0.)
        {
            return -PI/2.;
        }
        else if (y_ == 0)
        {
            // log_warn("If x and y are zero, the azimuth is not defined! Zero is returned!");
            return 0.;
        }
    }
    // log_fatal("never should be here; return zero.");
    return 0.;
}
