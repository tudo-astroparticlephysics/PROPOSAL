/*
 * Geometry.cxx
 *
 *  Created on: 05.06.2013
 *      Author: koehne
 */

#include "PROPOSAL/Geometry.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>

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

    if(inner_radius > radius_)
    {
        cerr<<"Warning: Inner radius is greater then radius (will be swaped)"<<endl;
        std::swap( inner_radius_ ,radius_ );
    }
    if(inner_radius == radius_)
    {
        cerr<<"Warning: Inner radius == radius (Volume is 0)"<<endl;
    }
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

    if(inner_radius > radius_)
    {
        cerr<<"Warning: Inner radius is greater then radius (will be swaped)"<<endl;
        std::swap( inner_radius_ ,radius_ );
    }
    if(inner_radius == radius_)
    {
        cerr<<"Warning: Inner radius == radius (Volume is 0)"<<endl;
    }

    z_      =   z;

    object_ =   "cylinder";
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Geometry::IsParticleInside(Particle* particle)
{
    bool is_inside  =   false;

    pair<double,double> dist =   DistanceToBorder(particle);

    if(dist.first > 0 && dist.second < 0)
    {
        is_inside   =   true;
    }
    return is_inside;


//    if( object_.compare("box")==0 )
//    {
//        double upper_x =   x0_ + 0.5*x_;
//        double lower_x =   x0_ - 0.5*x_;

//        double upper_y =   y0_ + 0.5*y_;
//        double lower_y =   y0_ - 0.5*y_;

//        double upper_z =   z0_ + 0.5*z_;
//        double lower_z =   z0_ - 0.5*z_;

//        //Figure out what happens if x/y/z == upper/lower...
//        if(    particle->GetX() < upper_x
//            && particle->GetX() > lower_x
//            && particle->GetY() < upper_y
//            && particle->GetY() > lower_y
//            && particle->GetZ() < upper_z
//            && particle->GetZ() > lower_z   )
//        {
//            is_inside   =   true;
//        }

//    }
//    else if( object_.compare("cylinder")==0 )
//    {
//        dist =    pow( (particle->GetX() -x0_) , 2.)
//                + pow( (particle->GetY() -y0_) , 2.);

//        double upper    =   z0_ + 0.5*z_;
//        double lower    =   z0_ - 0.5*z_;

//        //Figure out what happens if dist == radius/inner_radius z == upper/lower...
//        if(    particle->GetZ() < upper
//            && particle->GetZ() > lower
//            && dist             < radius_
//            && dist             > inner_radius_  )
//        {
//            is_inside   =   true;
//        }
//    }
//    else if( object_.compare("sphere")==0 )
//    {
//        dist =    pow( (particle->GetX() -x0_) , 2.)
//                + pow( (particle->GetY() -y0_) , 2.)
//                + pow( (particle->GetZ() -z0_) , 2.);

//        dist =  sqrt(dist);

//        //Figure out what happens if dist == radius/inner_radius...
//        if( dist < radius_ && dist > inner_radius_)
//        {
//            is_inside   =   true;
//        }

//    }
//    else
//    {
//        cerr<<"Error: In Geometry::IsParticleInside object name must be box/shpere/cylinder"<<endl;
//    }
//    return is_inside;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


pair<double,double> Geometry::DistanceToBorder(Particle* particle)
{
    //if( !IsParticleInside(particle) ) return 0;

    pair<double,double> distance;

    if( object_.compare("sphere")==0 )
    {
        distance    =   DistanceToBorderSphere(particle);

    }
    else if( object_.compare("box")==0 )
    {
        distance    =   DistanceToBorderBox(particle);

    }
//    else if( object_.compare("cylinder")==0 )
//    {
//        distance    =   DistanceToBorderCylinder(particle);

//    }
//    else
//    {
//        cerr<<"Warning: geometry type is not recognized! -1 is returned"<<endl;
//        return -1;
//    }

    return distance;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Geometry::Geometry()
    :x0_            ( 0. )
    ,y0_            ( 0. )
    ,z0_            ( 0. )
    ,inner_radius_  ( 0. )
    ,radius_        ( 0. )
    ,x_             ( 0. )
    ,y_             ( 0. )
    ,z_             ( 0. )
    ,object_        ( "sphere" )
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


pair<double,double> Geometry::DistanceToBorderSphere(Particle* particle)
{
    // Calculate intersection of particle trajectory and the sphere
    // sphere (x1 + x0)^2 + (x2 + y0)^2 + (x3 + z0)^2 = radius^2
    // straight line (particle trajectory) g = vec(x,y,z) + t * dir_vec( cosph *sinth, sinph *sinth , costh)
    // Insert and transform leads to C * t^2 + B * t + A = 0
    // length of direction vector =1 => C = 1
    // We are only interested in postive values of t
    // ( we want to find the intersection in direction of the particle trajectory)

    double A,B,t1,t2;
    double dir_vec_x = particle->GetCosPhi()*particle->GetSinTheta();
    double dir_vec_y = particle->GetSinPhi()*particle->GetSinTheta();
    double dir_vec_z = particle->GetCosTheta();

    pair<double,double> distance ;

    double determinant;

    A   =    pow( (particle->GetX() - x0_),2 )
           + pow( (particle->GetY() - y0_),2 )
           + pow( (particle->GetZ() - z0_),2 )
           - pow( radius_, 2 );

    B   =   2*(   (particle->GetX() - x0_)*dir_vec_x
                + (particle->GetY() - y0_)*dir_vec_y
                + (particle->GetZ() - z0_)*dir_vec_z );

    determinant =   pow(B/2 ,2) - A;

    if( determinant > 0) // determinant == 0 (boundery point) is ignored
    {
        t1  =   -1*B/2 + sqrt( determinant );
        t2  =   -1*B/2 - sqrt( determinant );

        // (-1/-1) sphere is behind particle or particle is on border but moving outside
        // ( dist_1 / dist_2 ) sphere is infront of the particle
        // ( dist_1 / -1 ) particle is inside the sphere or on border and moving inside
        if(t1 <= 0)
            distance.first  =   -1;

        else
            distance.first  =   t1;

        if(t2 <= 0)
            distance.second =   -1;

        else
            distance.second =   t2;

        // distance.first should be the smaller one
        if( distance.first < 0)
            std::swap(distance.first , distance.second);
        if( distance.first > 0 && distance.second > 0 )
        {
            if( distance.second < distance.first)
            {
                std::swap(distance.first , distance.second);
            }
        }


    }
    else    // particle trajectory does not have an intersection with the sphere
    {
        distance.first  =   -1;
        distance.second =   -1;
    }

    // No intersection so we don't have to check the distance to the inner sphere
    // if there is any
    if(distance.first < 0 && distance.second < 0)
        return distance;

    // This sqhere might be hollow and we have to check if the inner border is
    // reached before.
    // So we caluculate the intersection with the inner sphere.

    if(inner_radius_ > 0)
    {
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

            // Ok we have an intersection with the inner sphere

            // If distance.first and distance.second are positive this means
            // the sphere is infornt of the particle. So the first distance
            // ( intersection with the outer border) does not change
            // but the second distance has to be updated (intersection with the inner border)
            if(distance.first > 0 && distance.second > 0)
            {
                if( t1 > 0 )
                {
                    if( t1 < distance.second )
                        distance.second =   t1;

                }
                if( t2 > 0 )
                {
                    if( t2 < distance.second )
                        distance.second =   t2;

                }
            }
            else  // The particle is inside the outer sphere
            {
                //The particle is inside the inner sphere
                // which means outside the geometry
                if((t1 > 0 && t2 < 0) || (t2 > 0 && t1 < 0))
                {
                    distance.first  =   -1;
                }
                // Now we have to check if the particle is on the border of
                // the inner sphere
                if(t1 == 0 )
                {
                    // The particle is moving into the inner sphere
                    if( t2 > 0 )
                    {
                        distance.first  =   -1;
                    }
                    // if not we don't have to update distance.first
                }
                if(t2 == 0 )
                {
                    // The particle is moving into the inner sphere
                    if( t1 > 0 )
                    {
                        distance.first  =   -1;
                    }
                    // if not we don't have to update distance.first
                }
            }
        }
    }

    return distance;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


pair<double,double> Geometry::DistanceToBorderBox(Particle* particle)
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

    pair<double,double> distance;
    double t;
    double intersection_x;
    double intersection_y;
    double intersection_z;

    vector<double> dist;

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
                dist.push_back( t );
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
                dist.push_back( t );
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
                dist.push_back( t );
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
                dist.push_back( t );
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
                dist.push_back( t );
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
                dist.push_back( t );
            }

        }
    }

    if( dist.size() < 1 )   // No intersection with the box
    {
        distance.first  =   -1;
        distance.second =   -1;
    }
    else if( dist.size() == 1 )  // Particle is inside the box and we have one intersection in direction of the particle trajectory
    {
        distance.first  =   dist.at(0);
        distance.second =   -1;
    }
    else if( dist.size() ==2 )   // Particle is outside and the box is infront of the particle trajectory ( two intersections). Chose closest
    {
        distance.first  =   dist.at(0);
        distance.second =   dist.at(1);
        if( distance.second <   distance.first)
        {
            std::swap(distance.first , distance.second);
        }

    }
    else
    {
        cerr<<"In Geometry::DistanceToBorderCylinder(Particle*) this point should nerver be reached..."<<endl;
        cerr<<"(-1/-1) is returned"<<endl;

        distance.first  =   -1;
        distance.second =   -1;
    }

    return distance;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Geometry::DistanceToBorderCylinder(Particle* particle)
{
    // Calculate intersection of particle trajectory and the cylinder
    // cylinder barrel (x1 + x0)^2 + (x2 + y0)^2  = radius^2 [ z0_-0.5*z_ < particle->z <z0_ - 0.5*z_ ]
    // top/bottom surface:
    // E1: x3   =   z0_ + 0.5*z
    // E2: x3   =   z0_ - 0.5*z
    // straight line (particle trajectory) g = vec(x,y,z) + t * dir_vec( cosph *sinth, sinph *sinth , costh)
    // Insert and transform leads to C * t^2 + B * t + A = 0
    // length of direction vector =1 => C = 1
    // We are only interested in postive values of t
    // ( we want to find the intersection in direction of the particle trajectory)

    double A,B,t1,t2,t;
    double dir_vec_x = particle->GetCosPhi()*particle->GetSinTheta();
    double dir_vec_y = particle->GetSinPhi()*particle->GetSinTheta();
    double dir_vec_z = particle->GetCosTheta();

    double intersection_z;

    double distance =   0;

    if(!(dir_vec_x  == 0 && dir_vec_y == 0)) //Otherwise the particle trajectory is parallel to cylinder barrel
    {
        A   =    pow( (particle->GetX() - x0_),2 )
               + pow( (particle->GetY() - y0_),2 )
               - pow( radius_, 2 );

        B   =   2*(   (particle->GetX() - x0_)*dir_vec_x
                    + (particle->GetY() - y0_)*dir_vec_y );

        t1  =   -1*B/2 + sqrt( pow(B/2 ,2) - A );
        t2  =   -1*B/2 - sqrt( pow(B/2 ,2) - A );

        if(t1 > 0)
            t   =   t1;
        else
            t   =   t2;

        // t is the mupltiple of the direction vector we have to propagate
        // untill the border is reached.
        // But this cylinder might be hollow and we have to check if the inner border is
        // reached before.
        // So we caluculate the intersection with the inner sphere.

        if(inner_radius_ > 0)
        {
            double determinant;

            A   =    pow( (particle->GetX() - x0_),2 )
                   + pow( (particle->GetY() - y0_),2 )
                   - pow( radius_, 2 );

            B   =   2*(   (particle->GetX() - x0_)*dir_vec_x
                        + (particle->GetY() - y0_)*dir_vec_y );

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

        // Ok, now we know the multiple of the direction vector we have to propagte
        // untill the cylinder barrel is reached. But in z-direction this migth be outside
        // the cylinder border ( intersection with top or bottom surface )

        intersection_z  =   particle->GetZ() + t * dir_vec_z;

        // is inside the borders
        if( intersection_z > z0_ - 0.5*z_ &&
            intersection_z < z0_ - 0.5*z_    )
        {
            distance    =   t;     // Cause length of direction vector is 1 distance is t
            return distance;
        }

    }

    // When this point is reached ther must be an intersection with top or bottom surface

    //intersection with E1
    if( dir_vec_z != 0) // if dir_vec == 0 particle trajectory is parallel to E1 (Should not happen)
    {
        t   =  ( z0_ + 0.5*z_ - particle->GetZ() ) / dir_vec_z;

        if( t>0 ) // Interection is in particle trajectory direction
        {
            // We don't have to check for cylinder borders
            // cause if this intersection would not be inside the cylinder top or bottom surface
            // we found the intersection before when searching the intersection
            // with the cylinder barrel (the same is true for inner_radius > 0)

            distance    =   t;
            return  distance;

        }
    }

    //intersection with E2
    if( dir_vec_z != 0) // if dir_vec == 0 particle trajectory is parallel to E1 (Should not happen)
    {
        t   =  ( z0_ - 0.5*z_ - particle->GetZ() ) / dir_vec_z;

        if( t>0 ) // Interection is in particle trajectory direction
        {
            // We don't have to check for cylinder borders
            // cause if this intersection would not be inside the cylinder top or bottom surface
            // we found the intersection before when searching the intersection
            // with the cylinder barrel (the same is true for inner_radius > 0)

            distance    =   t;
            return  distance;

        }
    }

    cerr<<"In Geometry::DistanceToBorderCylinder(Particle*) this point should nerver be reached..."<<endl;
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
