#pragma once

#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <math.h>

namespace PROPOSAL{

class Vector3D
{
    public:
        // constructors
        Vector3D();
        Vector3D(const int x, const int y, const int z);
        Vector3D(const Vector3D& vector_3d);
        ~Vector3D(){};

        // operator functions and swap 
        Vector3D& operator=(const Vector3D& vector_3d);
        bool operator==(const Vector3D& vector_3d) const;
        bool operator!=(const Vector3D& vector_3d) const;
        void swap(Vector3D& vector_3d);

        //basic arithmetic
        Vector3D operator+ (const Vector3D& vec1, const Vector3D& vec2);
        Vector3D operator- (const Vector3D& vec1, const Vector3D& vec2);
        Vector3D operator* (const double factor1, const Vector3D& vec1);
        double scalar_product (const Vector3D& vec1, const Vector3D& vec2);
        Vector3D vector_product (const Vector3D& vec1, const Vector3D& vec2);
        double magnitude (const Vector3D& vec1, const Vector3D& vec2);

    private:
        double x_, y_, z_;
};

}

#endif