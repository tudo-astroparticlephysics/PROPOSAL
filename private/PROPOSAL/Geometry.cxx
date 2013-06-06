/*
 * Geometry.cxx
 *
 *  Created on: 05.06.2013
 *      Author: koehne
 */

#include "PROPOSAL/Geometry.h"
#include <iostream>
#include <cmath>

using namespace std;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------public member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Geometry::InitBox(double x0, double y0, double z0, double x, double y, double z)
{
    x0_     =   x0;
    y0_     =   y0;
    z0_     =   z0;

    x_      =   x;
    y_      =   y;
    z_      =   z;

    object_ =   "box";
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Geometry::InitSphere(double x0, double y0, double z0, double radius, double inner_radius)
{
    x0_     =   x0;
    y0_     =   y0;
    z0_     =   z0;

    radius_         =   radius;
    inner_radius_   =   inner_radius;

    object_ =   "sphere";
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Geometry::InitCylinder(double x0, double y0, double z0, double radius, double inner_radius, double z)
{
    x0_     =   x0;
    y0_     =   y0;
    z0_     =   z0;

    radius_         =   radius;
    inner_radius_   =   inner_radius;

    z_      =   z;

    object_ =   "cylinder";
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Geometry::IsParticleInside(Particle* particle)
{
    bool is_inside  =   false;

    double dist     =   0;

    if( object_.compare("box")==0 )
    {
        double upper_x =   x0_ + 0.5*x_;
        double lower_x =   x0_ - 0.5*x_;

        double upper_y =   y0_ + 0.5*y_;
        double lower_y =   y0_ - 0.5*y_;

        double upper_z =   z0_ + 0.5*z_;
        double lower_z =   z0_ - 0.5*z_;

        //Figure out what happens if x/y/z == upper/lower...
        if(    particle->GetX() < upper_x
            && particle->GetX() > lower_x
            && particle->GetY() < upper_y
            && particle->GetY() > lower_y
            && particle->GetZ() < upper_z
            && particle->GetZ() > lower_z   )
        {
            is_inside   =   true;
        }

    }
    else if( object_.compare("cylinder")==0 )
    {
        dist =    pow( (particle->GetX() -x0_) , 2.)
                + pow( (particle->GetY() -y0_) , 2.);

        double upper    =   z0_ + 0.5*z_;
        double lower    =   z0_ - 0.5*z_;

        //Figure out what happens if dist == radius/inner_radius z == upper/lower...
        if(    particle->GetZ() < upper
            && particle->GetZ() > lower
            && dist             < radius_
            && dist             > inner_radius_  )
        {
            is_inside   =   true;
        }
    }
    else if( object_.compare("sphere")==0 )
    {
        dist =    pow( (particle->GetX() -x0_) , 2.)
                + pow( (particle->GetY() -y0_) , 2.)
                + pow( (particle->GetZ() -z0_) , 2.);

        dist =  sqrt(dist);

        //Figure out what happens if dist == radius/inner_radius...
        if( dist < radius_ && dist > inner_radius_)
        {
            is_inside   =   true;
        }

    }
    else
    {
        cerr<<"Error: In Geometry::IsParticleInside object name must be box/shpere/cylinder"<<endl;
    }
    return is_inside;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Geometry::DistanceToBorder(Particle* particle)
{
    if( !IsParticleInside(particle) ) return 0;

    double distance =   0;

    if( object_.compare("sphere")==0 )
    {
        distance    =   DistanceToBorderSphere(particle);

    }


    return distance;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Geometry::Geometry()
    :x0_            ( 0 )
    ,y0_            ( 0 )
    ,z0_            ( 0 )
    ,inner_radius_  ( 0 )
    ,radius_        ( 0 )
    ,x_             ( 0 )
    ,y_             ( 0 )
    ,z_             ( 0 )
    ,object_        ( "sqhere" )
{

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Geometry::Geometry(const Geometry &geometry)
    :x0_            ( geometry.x0_ )
    ,y0_            ( geometry.y0_ )
    ,z0_            ( geometry.z0_ )
    ,inner_radius_  ( geometry.inner_radius_ )
    ,radius_        ( geometry.radius_ )
    ,x_             ( geometry.x_ )
    ,y_             ( geometry.y_ )
    ,z_             ( geometry.z_ )
    ,object_        ( geometry.object_ )
{

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Geometry& Geometry::operator=(const Geometry &geometry)
{
    if (this != &geometry)
    {
      Geometry tmp(geometry);
      swap(tmp);
    }
    return *this;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Geometry::operator==(const Geometry &geometry) const
{
    if( x0_            != geometry.x0_ )            return false;
    if( y0_            != geometry.y0_ )            return false;
    if( z0_            != geometry.z0_ )            return false;
    if( inner_radius_  != geometry.inner_radius_ )  return false;
    if( radius_        != geometry.radius_ )        return false;
    if( x_             != geometry.x_ )             return false;
    if( y_             != geometry.y_ )             return false;
    if( z_             != geometry.z_ )             return false;

    if( object_.compare(geometry.object_)!= 0 )     return false;

    //else
    return true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Geometry::operator!=(const Geometry &geometry) const
{
    return !(*this == geometry);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Geometry::swap(Geometry &geometry)
{
    using std::swap;

    swap( x0_            , geometry.x0_ );
    swap( y0_            , geometry.y0_ );
    swap( z0_            , geometry.z0_ );
    swap( inner_radius_  , geometry.inner_radius_ );
    swap( radius_        , geometry.radius_ );
    swap( x_             , geometry.x_ );
    swap( y_             , geometry.y_ );
    swap( z_             , geometry.z_ );

    object_.swap( geometry.object_ );
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------private member functions---------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Geometry::DistanceToBorderSphere(Particle* particle)
{
    // Calculate intersection of particle trajectory and the sphere
    // sphere (x1 + x0)^2 + (x2 + y0)^2 + (x3 + z0)^2 = radius^2
    // straight line (particle trajectory) g = vec(x,y,z) + t * dir_vec( cosph *sinth, sinph *sinth , costh)
    // Insert and transform leads to C * t^2 + B * t + A = 0
    // length of direction vector =1 => C = 1
    // We are only interested in postive values of t
    // ( we want to find the intersection in direction of the particle trajectory)

    double A,B,t1,t2,t;
    double dir_vec_x = particle->GetCosPhi()*particle->GetSinTheta();
    double dir_vec_y = particle->GetSinPhi()*particle->GetSinTheta();
    double dir_vec_z = particle->GetCosTheta();

    double distance =   0;


    A   =    pow( (particle->GetX() - x0_),2 )
           + pow( (particle->GetY() - y0_),2 )
           + pow( (particle->GetZ() - z0_),2 )
           - pow( radius_, 2 );

    B   =   2*(   (particle->GetX() - x0_)*dir_vec_x
                + (particle->GetY() - y0_)*dir_vec_y
                + (particle->GetZ() - z0_)*dir_vec_z );

    t1  =   -1*B/2 + sqrt( pow(B/2 ,2) - A );
    t2  =   -1*B/2 - sqrt( pow(B/2 ,2) - A );

    if(t1 > 0)
        t   =   t1;
    else
        t   =   t2;

    // t is the mupltiple of the direction vector we have to propagate untill the
    // border is reached.
    // But this sqhere might be hollow and we have to check if the inner border is
    // reached before.
    // So we caluculate the intersection with the inner sphere.

    if(inner_radius_ > 0)
    {
        double determinant;

        A   =    pow( (particle->GetX() - x0_),2 )
               + pow( (particle->GetY() - y0_),2 )
               + pow( (particle->GetZ() - z0_),2 )
               - pow( inner_radius_, 2 );

        B   =   2*(   (particle->GetX() - x0_)*dir_vec_x
                    + (particle->GetY() - y0_)*dir_vec_y
                    + (particle->GetZ() - z0_)*dir_vec_z );

        determinant = pow(B/2 ,2) - A;

        if( determinant > 0) // determinant == 0 (boundery point) is ignored
        {
            t1  =   -1*B/2 + sqrt( determinant );
            t2  =   -1*B/2 - sqrt( determinant );

            // Ok we have a intersection with the inner sphere
            // Now we have to find the closest
            if( t1 != 0 && t1 < t)
            {
                t   =   t1;
            }
            if( t2 != 0 && t2 < t)
            {
                t   =   t2;
            }
        }
    }

    //Cause length of direction vector is 1 distance is t

    distance    =   t;

    return distance;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Geometry::DistanceToBorderBox(Particle* particle)
{
    // Calculate intersection of particle trajectory and the box
    // Surface of the box is defined by six planes:
    // E1: x1   =   x0_ + 0.5*x
    // E2: x1   =   x0_ - 0.5*x
    // E3: x2   =   y0_ + 0.5*y
    // E4: x2   =   y0_ - 0.5*y
    // E5: x3   =   z0_ + 0.5*z
    // E6: x3   =   z0_ - 0.5*z
    // straight line (particle trajectory) g = vec(x,y,z) + t * dir_vec( cosph *sinth, sinph *sinth , costh)
    // We are only interested in postive values of t
    // ( we want to find the intersection in direction of the particle trajectory)

    double dir_vec_x = particle->GetCosPhi()*particle->GetSinTheta();
    double dir_vec_y = particle->GetSinPhi()*particle->GetSinTheta();
    double dir_vec_z = particle->GetCosTheta();

    double distance =   0;
    double t;
    double intersection_x;
    double intersection_y;
    double intersection_z;

    //intersection with E1
    if( dir_vec_x != 0) // if dir_vec == 0 particle trajectory is parallel to E1
    {
        t   =  ( x0_ + 0.5*x_ - particle->GetX() ) / dir_vec_x;

        if( t>0 ) // Interection is in particle trajectory direction
        {
            //Check if intersection is inside the box borders
            intersection_y  =   particle->GetY() + t * dir_vec_y;
            intersection_z  =   particle->GetZ() + t * dir_vec_z;
            if( intersection_y >= (y0_ - 0.5*y_) &&
                intersection_y <= (y0_ + 0.5*y_) &&
                intersection_z >= (z0_ - 0.5*z_) &&
                intersection_z <= (z0_ + 0.5*z_)    )
            {
                distance    =   t;
                return  distance;
            }

        }
    }

    //intersection with E2
    if( dir_vec_x != 0) // if dir_vec == 0 particle trajectory is parallel to E2
    {
        t   =  ( x0_ - 0.5*x_ - particle->GetX() ) / dir_vec_x;

        if( t>0 ) // Interection is in particle trajectory direction
        {
            //Check if intersection is inside the box borders
            intersection_y  =   particle->GetY() + t * dir_vec_y;
            intersection_z  =   particle->GetZ() + t * dir_vec_z;
            if( intersection_y >= (y0_ - 0.5*y_) &&
                intersection_y <= (y0_ + 0.5*y_) &&
                intersection_z >= (z0_ - 0.5*z_) &&
                intersection_z <= (z0_ + 0.5*z_)    )
            {
                distance    =   t;
                return  distance;
            }

        }
    }

    //intersection with E3
    if( dir_vec_y != 0) // if dir_vec == 0 particle trajectory is parallel to E3
    {
        t   =  ( y0_ + 0.5*y_ - particle->GetY() ) / dir_vec_y;

        if( t>0 ) // Interection is in particle trajectory direction
        {
            //Check if intersection is inside the box borders
            intersection_x  =   particle->GetX() + t * dir_vec_x;
            intersection_z  =   particle->GetZ() + t * dir_vec_z;
            if( intersection_x >= (x0_ - 0.5*x_) &&
                intersection_x <= (x0_ + 0.5*x_) &&
                intersection_z >= (z0_ - 0.5*z_) &&
                intersection_z <= (z0_ + 0.5*z_)    )
            {
                distance    =   t;
                return  distance;
            }

        }
    }

    //intersection with E4
    if( dir_vec_y != 0) // if dir_vec == 0 particle trajectory is parallel to E4
    {
        t   =  ( y0_ - 0.5*y_ - particle->GetY() ) / dir_vec_y;

        if( t>0 ) // Interection is in particle trajectory direction
        {
            //Check if intersection is inside the box borders
            intersection_x  =   particle->GetX() + t * dir_vec_x;
            intersection_z  =   particle->GetZ() + t * dir_vec_z;
            if( intersection_x >= (x0_ - 0.5*x_) &&
                intersection_x <= (x0_ + 0.5*x_) &&
                intersection_z >= (z0_ - 0.5*z_) &&
                intersection_z <= (z0_ + 0.5*z_)    )
            {
                distance    =   t;
                return  distance;
            }

        }
    }

    //intersection with E5
    if( dir_vec_z != 0) // if dir_vec == 0 particle trajectory is parallel to E5
    {
        t   =  ( z0_ + 0.5*z_ - particle->GetZ() ) / dir_vec_z;

        if( t>0 ) // Interection is in particle trajectory direction
        {
            //Check if intersection is inside the box borders
            intersection_x  =   particle->GetX() + t * dir_vec_x;
            intersection_y  =   particle->GetY() + t * dir_vec_y;
            if( intersection_x >= (x0_ - 0.5*x_) &&
                intersection_x <= (x0_ + 0.5*x_) &&
                intersection_y >= (y0_ - 0.5*y_) &&
                intersection_y <= (y0_ + 0.5*y_)    )
            {
                distance    =   t;
                return  distance;
            }

        }
    }

    //intersection with E6
    if( dir_vec_z != 0) // if dir_vec == 0 particle trajectory is parallel to E6
    {
        t   =  ( z0_ - 0.5*z_ - particle->GetZ() ) / dir_vec_z;

        if( t>0 ) // Interection is in particle trajectory direction
        {
            //Check if intersection is inside the box borders
            intersection_x  =   particle->GetX() + t * dir_vec_x;
            intersection_y  =   particle->GetY() + t * dir_vec_y;
            if( intersection_x >= (x0_ - 0.5*x_) &&
                intersection_x <= (x0_ + 0.5*x_) &&
                intersection_y >= (y0_ - 0.5*y_) &&
                intersection_y <= (y0_ + 0.5*y_)    )
            {
                distance    =   t;
                return  distance;
            }

        }
    }

    cerr<<"In Geometry::DistanceToBorderBox(Particle*) this point should nerver be reached..."<<endl;
    cerr<<"-1 is returned"<<endl;

    return -1;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Setter-------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Geometry::SetX0(double x0) {
    x0_ =   x0;
}

void Geometry::SetY0(double y0) {
    y0_ =   y0;
}

void Geometry::SetZ0(double z0) {
    z0_ =   z0;
}

void Geometry::SetX(double x) {
    x_ =   x;
}

void Geometry::SetY(double y) {
    y_ =   y;
}

void Geometry::SetZ(double z) {
    z_ =   z;
}

void Geometry::SetInnerRadius(double inner_radius_) {
    inner_radius_ =   inner_radius_;
}

void Geometry::SetRadius(double radius) {
    radius_ =   radius;
}

void Geometry::SetObject(string object) {
    object_ =   object;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Destructor---------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Geometry::~Geometry(){}
