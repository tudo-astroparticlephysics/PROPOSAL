

#include "PROPOSAL/Vector3D.h"

using namespace PROPOSAL;


//----------------------------------------------------------------------//
//------------------------- Constructors -------------------------------//
//----------------------------------------------------------------------//

// default constructor
Vector3D::Vector3D()
    : x_ (0)
    , y_ (0)
    , z_ (0)
{
}


Vector3D::Vector3D(const double x, const double y, const double z)
    : x_ (x)
    , y_ (y)
    , z_ (z)
{
}

// copy constructor
Vector3D::Vector3D(const Vector3D& vector_3d)
    : x_ (vector_3d.x_)
    , y_ (vector_3d.y_)
    , z_ (vector_3d.z_)
{
}

// destructor
Vector3D::~Vector3D(){}

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
    if(x_ != vector_3d.x_) return false;
    if(y_ != vector_3d.y_) return false;
    if(z_ != vector_3d.z_) return false;

    return true;
}

bool Vector3D::operator!=(const Vector3D& vector_3d) const
{
    return !(*this == vector_3d)
}

void Vector3D::swap(Vector3D& vector_3d)
{
    using std::swap;

    swap(x_, vector_3d.x_);
    swap(y_, vector_3d.y_);
    swap(z_, vector_3d.z_);
}

//----------------------------------------------------------------------//
//-----------------------operator basic arithmetic ---------------------//
//----------------------------------------------------------------------//

Vector3D Vector3D::operator+(const Vector3D& vec1, const Vector3D& vec2)
{
    Vector3D vector_sum;
    vector_sum.x_ = vec1.x_ + vec2.x_;
    vector_sum.y_ = vec1.y_ + vec2.y_;
    vector_sum.z_ = vec1.z_ + vec2.z_;
    return vector_sum;
}

Vector3D Vector3D::operator-(const Vector3D& vec1, const Vector3D& vec2)
{
    Vector3D vector_diff;
    vector_diff.x_ = vec1.x_ - vec2.x_;
    vector_diff.y_ = vec1.y_ - vec2.y_;
    vector_diff.z_ = vec1.z_ - vec2.z_;
    return vector_diff;
}

Vector3D Vector3D::operator*(const double factor1, const Vector3D& vec1)
{
    Vector3D product;
    product.x_ = factor1*vec1.x_;
    product.y_ = factor1*vec1.y_;
    product.z_ = factor1*vec1.z_;
    return vector_product;
}

double Vector3D::scalar_product(const Vector3D& vec1, const Vector3D& vec2)
{
    return vec1.x_*vec2.x_ + vec1.y_*vec2.y_ + vec1.z_*vec2.z_;
}

Vector3D Vector3D::vector_product(const Vector3D& vec1, const Vector3D& vec2)
{
    Vector3D product;
    product.x_ = vec1.y_*vec2.z_ - vec1.z_*vec2.y_;
    product.y_ = vec1.z_*vec2.x_ - vec1.x_*vec2.z_;
    product.z_ = vec1.x_*vec2.y_ - vec1.y_*vec2.x_;
    return product;
}

double Vector3D::magnitude(const Vector3D& vector_3d)
{
    double aux;
    aux = vector_3d.x_*vector_3d.x_
        + vector_3d.y_*vector_3d.y_
        + vector_3d.z_*vector_3d.z_;
    return aux;
}

