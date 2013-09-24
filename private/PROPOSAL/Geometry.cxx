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
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Output.h"

using namespace std;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------public member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Geometry::InitBox(double x0, double y0, double z0, double x, double y, double z)
{
    x0_     =   100*x0;
    y0_     =   100*y0;
    z0_     =   100*z0;

    x_      =   100*x;
    y_      =   100*y;
    z_      =   100*z;

    object_ =   "box";
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Geometry::InitSphere(double x0, double y0, double z0, double radius, double inner_radius)
{
    x0_     =   100*x0;
    y0_     =   100*y0;
    z0_     =   100*z0;

    radius_         =   100*radius;
    inner_radius_   =   100*inner_radius;

    if(inner_radius > radius_)
    {
        log_error("Inner radius %f is greater then radius %f (will be swaped)",inner_radius_ ,radius_);
        std::swap( inner_radius_ ,radius_ );
    }
    if(inner_radius == radius_)
    {
        log_error("Warning: Inner radius %f == radius %f (Volume is 0)",inner_radius_,radius_);
    }
    object_ =   "sphere";
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Geometry::InitCylinder(double x0, double y0, double z0, double radius, double inner_radius, double z)
{
    x0_     =   100*x0;
    y0_     =   100*y0;
    z0_     =   100*z0;

    radius_         =   100*radius;
    inner_radius_   =   100*inner_radius;

    if(inner_radius_ > radius_)
    {
        log_error("Inner radius %f is greater then radius %f (will be swaped)",inner_radius_ ,radius_);
        std::swap( inner_radius_ ,radius_ );
    }
    if(inner_radius_ == radius_)
    {
        log_error("Warning: Inner radius %f == radius %f (Volume is 0)",inner_radius_,radius_);
    }

    z_      =   100*z;

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
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Geometry::IsParticleInfront(Particle* particle)
{
    bool is_infront  =   false;

    pair<double,double> dist =   DistanceToBorder(particle);

    if(dist.first > 0 && dist.second > 0)
    {
        is_infront   =   true;
    }
    return is_infront;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Geometry::IsParticleBehind(Particle* particle)
{
    bool is_behind  =   false;

    pair<double,double> dist =   DistanceToBorder(particle);

    if(dist.first < 0 && dist.second < 0)
    {
        is_behind   =   true;
    }
    return is_behind;
}


//----------------------------------------------------------------------------//
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
    else if( object_.compare("cylinder")==0 )
    {
        distance    =   DistanceToBorderCylinder(particle);

    }
    else
    {
        log_error("geometry type is not recognized! (-1,-1) is returned");
        distance.first  =   -1;
        distance.second =   -1;
    }

    return distance;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


double Geometry::DistanceToClosestApproach(Particle* particle)
{
    double dir_vec_x = particle->GetCosPhi()*particle->GetSinTheta();
    double dir_vec_y = particle->GetSinPhi()*particle->GetSinTheta();
    double dir_vec_z = particle->GetCosTheta();

    double distance  =   (x0_ - particle->GetX()) * dir_vec_x +
                         (y0_ - particle->GetY()) * dir_vec_y +
                         (z0_ - particle->GetZ()) * dir_vec_z;
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
    ,hirarchy_      (0)
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
    ,hirarchy_      ( geometry.hirarchy_)
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

    if( hirarchy_      != geometry.hirarchy_  )     return false;

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


ostream& operator<<(ostream& os, Geometry const& geometry)
{
    os<<"--------Geometry( "<<&geometry<<" )--------"<<endl;
    os<<"\t"<<geometry.object_<<endl;
    os<<"\tOrigin (x,y,z):\t"<<geometry.x0_<<"\t"<<geometry.y0_<<"\t"<<geometry.z0_<<endl;
    os<<"\tHirarchy:\t"<<geometry.hirarchy_ << endl;
    if(geometry.object_.compare("sphere")==0)
    {
        os<<"\tRadius: "<<geometry.radius_<<"\tInner radius: "<<geometry.inner_radius_<<endl;
    }
    if(geometry.object_.compare("cylinder")==0)
    {
        os<<"\tRadius: "<<geometry.radius_<<"\tInnner radius: "<<geometry.inner_radius_<<" Height: "<<geometry.z_<<endl;
    }
    if(geometry.object_.compare("box")==0)
    {
        os<<"\tWidth_x: "<<geometry.x_<<"\tWidth_y "<<geometry.y_<<"\tHeight: "<<geometry.z_<<endl;
    }
    os<<"------------------------------------";
    return os;
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
    swap( hirarchy_      , geometry.hirarchy_ );

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

        //Computer precision controll
        if(t1 > 0 && t1 < GEOMETRY_PRECISION)
            t1  =   0;
        if(t2 > 0 && t2 < GEOMETRY_PRECISION)
            t2  =   0;

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

            //Computer precision controll
            if(t1 > 0 && t1 < GEOMETRY_PRECISION)
                t1  =   0;
            if(t2 > 0 && t2 < GEOMETRY_PRECISION)
                t2  =   0;

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
                //The inner cylinder is infront of the particle trajectory
                //distance.first has to be updated
                if(t1 > 0 && t2 > 0)
                {
                    if(t1<t2)
                        distance.first  =   t1;
                    else
                        distance.first  =   t2;
                }
                // The particle is inside the inner sphere
                // this means distance.second becomes distanc.first
                // and distance.first beomces distance to intersection with
                // the inner sphere in direction of the particle trajectory
                if((t1 > 0 && t2 < 0) || (t2 > 0 && t1 < 0))
                {
                    std::swap(distance.first,distance.second);
                    if( t1 > 0 )
                        distance.first  =   t1;
                    else
                        distance.first  =   t2;
                }
                // Now we have to check if the particle is on the border of
                // the inner sphere
                if(t1 == 0 )
                {
                    // The particle is moving into the inner sphere
                    if( t2 > 0 )
                    {
                        std::swap(distance.first,distance.second);
                        distance.first  =   t2;
                    }
                    // if not we don't have to update distance.first
                }
                if(t2 == 0 )
                {
                    // The particle is moving into the inner sphere
                    if( t1 > 0 )
                    {
                        std::swap(distance.first,distance.second);
                        distance.first  =   t1;
                    }
                    // if not we don't have to update distance.first
                }
            }
        }
    }
    // Make a computer precision controll!
    // This is necessary cause due to numerical effects it meight be happen
    // that a particle which is located on a gemoetry border is treated as inside
    // or outside

    if( distance.first < GEOMETRY_PRECISION )
        distance.first  =   -1;
    if( distance.second < GEOMETRY_PRECISION )
        distance.second  =   -1;
    if( distance.first < 0 )
        std::swap(distance.first ,distance.second);

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

        //Computer precision controll
        if(t > 0 && t < GEOMETRY_PRECISION)
            t  =   0;

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

        //Computer precision controll
        if(t > 0 && t < GEOMETRY_PRECISION)
            t  =   0;

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

        //Computer precision controll
        if(t > 0 && t < GEOMETRY_PRECISION)
            t  =   0;

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

        //Computer precision controll
        if(t > 0 && t < GEOMETRY_PRECISION)
            t  =   0;

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

        //Computer precision controll
        if(t > 0 && t < GEOMETRY_PRECISION)
            t  =   0;

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

        //Computer precision controll
        if(t > 0 && t < GEOMETRY_PRECISION)
            t  =   0;

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
    else if( dist.size() ==2 )   // Particle is outside and the box is infront of the particle trajectory ( two intersections).
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
        log_error("This point should nerver be reached... (-1/-1) is returned");

        distance.first  =   -1;
        distance.second =   -1;
    }

    // Make a computer precision controll!
    // This is necessary cause due to numerical effects it meight be happen
    // that a particle which is located on a gemoetry border is treated as inside
    // or outside
    if( distance.first < GEOMETRY_PRECISION )
        distance.first  =   -1;
    if( distance.second < GEOMETRY_PRECISION )
        distance.second =   -1;
    if( distance.first < 0 )
        std::swap(distance.first ,distance.second);

    return distance;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


pair<double,double> Geometry::DistanceToBorderCylinder(Particle* particle)
{
    // Calculate intersection of particle trajectory and the cylinder
    // cylinder barrel (x1 + x0)^2 + (x2 + y0)^2  = radius^2 [ z0_-0.5*z_ < particle->z <z0_ - 0.5*z_ ]
    // top/bottom surface:
    // E1: x3   =   z0_ + 0.5*z
    // E2: x3   =   z0_ - 0.5*z
    // straight line (particle trajectory) g = vec(x,y,z) + t * dir_vec( cosph *sinth, sinph *sinth , costh)
    // Insert and transform leads to C * t^2 + B * t + A = 0
    // We are only interested in postive values of t
    // ( we want to find the intersection in direction of the particle trajectory)

    // (-1/-1) cylinder is behind particle or particle is on border but moving outside
    // ( dist_1 / dist_2 ) cylinder is infront of the particle
    // ( dist_1 / -1 ) particle is inside the cylinder or on border and moving inside

    double A,B,C,t1,t2,t;
    double dir_vec_x = particle->GetCosPhi()*particle->GetSinTheta();
    double dir_vec_y = particle->GetSinPhi()*particle->GetSinTheta();
    double dir_vec_z = particle->GetCosTheta();

    double determinant;

    double intersection_x;
    double intersection_y;
    double intersection_z;

    vector<double> dist;

    pair<double,double> distance;

    if(!(dir_vec_x  == 0 && dir_vec_y == 0)) //Otherwise the particle trajectory is parallel to cylinder barrel
    {

        A   =    pow( (particle->GetX() - x0_),2 )
               + pow( (particle->GetY() - y0_),2 )
               - pow( radius_, 2 );

        B   =   2*(   (particle->GetX() - x0_)*dir_vec_x
                    + (particle->GetY() - y0_)*dir_vec_y );

        C   =   dir_vec_x*dir_vec_x + dir_vec_y*dir_vec_y;

        B   /=  C;
        A   /=  C;


        determinant =   pow(B/2 ,2) - A;

        if( determinant > 0) // determinant == 0 (boundery point) is ignored
        {
            t1  =   -1*B/2 + sqrt( determinant );
            t2  =   -1*B/2 - sqrt( determinant );

            //Computer precision controll
            if(t1 > 0 && t1 < GEOMETRY_PRECISION)
                t1  =   0;
            if(t2 > 0 && t2 < GEOMETRY_PRECISION)
                t2  =   0;


            if(t1 > 0)
            {
                intersection_z  =   particle->GetZ() + t1 * dir_vec_z;
                // is inside the borders
                if( intersection_z > z0_ - 0.5*z_ &&
                    intersection_z < z0_ + 0.5*z_    )
                {
                    dist.push_back(t1);
                }

            }

            if(t2 > 0)
            {
                intersection_z  =   particle->GetZ() + t2 * dir_vec_z;
                // is inside the borders
                if( intersection_z > z0_ - 0.5*z_ &&
                    intersection_z < z0_ + 0.5*z_    )
                {
                    dist.push_back(t2);
                }
            }
        }
    }

    // if we have found already to intersections we don't have to check for intersections
    // with top or bottom surface
    if( dist.size() < 2 )
    {
        //intersection with E1
        if( dir_vec_z != 0) // if dir_vec == 0 particle trajectory is parallel to E1 (Should not happen)
        {
            t   =  ( z0_ + 0.5*z_ - particle->GetZ() ) / dir_vec_z;
            //Computer precision controll
            if(t > 0 && t < GEOMETRY_PRECISION)
                t  =   0;

            if( t>0 ) // Interection is in particle trajectory direction
            {
                intersection_x  =   particle->GetX() + t * dir_vec_x;
                intersection_y  =   particle->GetY() + t * dir_vec_y;

                if( sqrt( pow( (intersection_x -x0_) ,2) + pow( (intersection_y - y0_) ,2 ) ) <= radius_ &&
                    sqrt( pow( (intersection_x -x0_) ,2) + pow( (intersection_y - y0_) ,2 ) ) >= inner_radius_ )
                {
                    dist.push_back( t );
                }
            }
        }

        //intersection with E2
        if( dir_vec_z != 0) // if dir_vec == 0 particle trajectory is parallel to E2 (Should not happen)
        {
            t   =  ( z0_ - 0.5*z_ - particle->GetZ() ) / dir_vec_z;

            //Computer precision controll
            if(t > 0 && t < GEOMETRY_PRECISION)
                t  =   0;

            if( t>0 ) // Interection is in particle trajectory direction
            {
                intersection_x  =   particle->GetX() + t * dir_vec_x;
                intersection_y  =   particle->GetY() + t * dir_vec_y;

                if(sqrt( pow( (intersection_x - x0_) ,2) + pow( (intersection_y - y0_) ,2 ) ) <= radius_ &&
                   sqrt( pow( (intersection_x - x0_) ,2) + pow( (intersection_y - y0_) ,2 ) ) >= inner_radius_ )
                {
                    dist.push_back( t );
                }
            }
        }
    }
    // No intersection with the outer cylinder
    if(dist.size() < 1)
    {
        distance.first  =   -1;
        distance.second =   -1;
    //    return distance;
    }
    else if(dist.size() == 1)   //particle is inside the cylinder
    {
        distance.first =    dist.at(0);
        distance.second=    -1;

    }
    else if(dist.size() == 2)   //cylinder is infront of the particle
    {
        distance.first =    dist.at(0);
        distance.second=    dist.at(1);

        if(distance.second < distance.first)
        {
            std::swap(distance.first , distance.second);
        }

    }
    else
    {
        log_error("This point should never be reached");
    }
    // This cylinder might be hollow and we have to check if the inner border is
    // reached before.
    // So we caluculate the intersection with the inner cylinder.

    if(inner_radius_ > 0)
    {
        if(!(dir_vec_x  == 0 && dir_vec_y == 0))
        {

            A   =    pow( (particle->GetX() - x0_),2 )
                   + pow( (particle->GetY() - y0_),2 )
                   - pow( inner_radius_, 2 );

            B   =   2*(   (particle->GetX() - x0_)*dir_vec_x
                        + (particle->GetY() - y0_)*dir_vec_y );

            C   =   dir_vec_x*dir_vec_x + dir_vec_y*dir_vec_y;

            B   /=  C;
            A   /=  C;

            determinant = pow(B/2 ,2) - A;

            if( determinant > 0) // determinant == 0 (boundery point) is ignored
            {
                t1  =   -1*B/2 + sqrt( determinant );
                t2  =   -1*B/2 - sqrt( determinant );

                //Computer precision controll
                if(t1 > 0 && t1 < GEOMETRY_PRECISION)
                    t1  =   0;
                if(t2 > 0 && t2 < GEOMETRY_PRECISION)
                    t2  =   0;

                // Ok we have a intersection with the inner cylinder

                // If distance.first and distance.second are positive this means
                // the cylinder is infornt of the particle. So the first distance
                // ( intersection with the outer border) does not change
                // but the second distance has to be updated (intersection with the inner border)
                if(distance.first > 0 && distance.second > 0)
                {
                    if(t1 > 0)
                    {
                        intersection_z  =   particle->GetZ() + t1 * dir_vec_z;
                        // is inside the borders
                        if( intersection_z > z0_ - 0.5*z_ &&
                            intersection_z < z0_ + 0.5*z_    )
                        {
                            if( t1 < distance.second )
                                distance.second =   t1;
                        }

                    }

                    if(t2 > 0)
                    {
                        intersection_z  =   particle->GetZ() + t2 * dir_vec_z;
                        // is inside the borders
                        if( intersection_z > z0_ - 0.5*z_ &&
                            intersection_z < z0_ + 0.5*z_    )
                        {
                            if( t2 < distance.second )
                                distance.second =   t2;
                        }

                    }
                }
                // The particle trajectory hits the inner cylinder first
                // (particle flys through the hole and hits inner cylinder barrel first)
                //  ___  ^     ___
                // |   |  \   |   |
                // |   |   \  |   |
                // |   |    \ |   |
                // |   |     \|   |
                // |   |      *   |
                // |   |      |\  |
                // |   |      | x |
                // |   |      |   |
                // |___|      |___|
                else if(distance.first < 0 && distance.second < 0)
                {
                    if(t1 > 0)
                    {
                        intersection_z  =   particle->GetZ() + t1 * dir_vec_z;
                        // is inside the borders
                        if( intersection_z > z0_ - 0.5*z_ &&
                            intersection_z < z0_ + 0.5*z_    )
                        {
                                distance.first =   t1;
                        }

                    }
                    if(t2 > 0)
                    {
                        intersection_z  =   particle->GetZ() + t2 * dir_vec_z;
                        // is inside the borders
                        if( intersection_z > z0_ - 0.5*z_ &&
                            intersection_z < z0_ + 0.5*z_    )
                        {
                                distance.first =   t2;
                        }

                    }
                }
                // The particle is inside the outer cylinder
                // or particle hits bottom or top surface and then hits
                // the inner cylinder barrel
                else
                {
                    //The inner cylinder is infront of the particle trajectory
                    //distance.first has to be updated
                    if(t1 > 0 && t2 > 0)
                    {
                        //  _____        _____
                        // |     |      |     |
                        // |     |      |     |
                        // |     |      |     |
                        // |     |      |     |
                        // |     |      |     |
                        // | x---*----->|     |
                        // |     |      |     |
                        // |     |      |     |
                        // |_____|      |_____|
                        //
                        if( particle->GetZ() >= z0_ - 0.5*z_ &&
                            particle->GetZ() <= z0_ + 0.5*z_ &&
                            sqrt( pow( (particle->GetX() - x0_) ,2) + pow( (particle->GetY() - y0_) ,2 ) ) <= radius_ + GEOMETRY_PRECISION &&
                            sqrt( pow( (particle->GetX() - x0_) ,2) + pow( (particle->GetY() - y0_) ,2 ) ) >= inner_radius_ -GEOMETRY_PRECISION  )
                        {
                            if( t1 < distance.first )
                            {
                                intersection_z  =   particle->GetZ() + t1 * dir_vec_z;
                                // is inside the borders
                                if( intersection_z > z0_ - 0.5*z_ &&
                                    intersection_z < z0_ + 0.5*z_    )
                                {
                                    // This case means particle is inside in hits
                                    // inner cylinder first
                                    distance.first  =   t1;
                                }
                            }
                            if( t2 < distance.first )
                            {
                                intersection_z  =   particle->GetZ() + t2 * dir_vec_z;
                                // is inside the borders
                                if( intersection_z > z0_ - 0.5*z_ &&
                                    intersection_z < z0_ + 0.5*z_    )
                                {
                                    // This case means particle is inside in hits
                                    // inner cylinder first
                                    distance.first  =   t2;
                                }
                            }
                        }
                        //               ^
                        //  _____       / _____
                        // |     |     / |     |
                        // |     |    /  |     |
                        // |     |   /   |     |
                        // |     |  /    |     |
                        // |     | /     |     |
                        // |     |/      |     |
                        // |     *       |     |
                        // |    /|       |     |
                        // |___/_|       |_____|
                        //    *
                        //   x
                        else
                        {
                            intersection_z  =   particle->GetZ() + t1 * dir_vec_z;
                            // is inside the borders
                            if( intersection_z > z0_ - 0.5*z_ &&
                                intersection_z < z0_ + 0.5*z_    )
                            {
                                // This case means particle is inside in hits
                                // inner cylinder first
                                distance.second  =   t1;
                            }

                            if(distance.second < 0)
                            {
                                intersection_z  =   particle->GetZ() + t2 * dir_vec_z;
                                // is inside the borders
                                if( intersection_z > z0_ - 0.5*z_ &&
                                    intersection_z < z0_ + 0.5*z_    )
                                {
                                    // This case means particle is inside in hits
                                    // inner cylinder first
                                    distance.second  =   t2;
                                }
                            }
                        }
                    }
                    // The particle is inside the inner cylinder
                    // this means distance.second becomes distanc.first
                    // and distance.first beomces distance to intersection with
                    // the inner cylinder in direction of the particle trajectory
                    //  _____        _____
                    // |     |      |     |
                    // |     |      |     |
                    // |     |      |     |
                    // |     |      |     |
                    // |     |      |     |
                    // |     |  x---*-----*------->
                    // |     |      |     |
                    // |     |      |     |
                    // |_____|      |_____|
                    //
                    if((t1 > 0 && t2 < 0) || (t2 > 0 && t1 < 0))
                    {
                        if( t1 > 0 )
                        {
                            intersection_z  =   particle->GetZ() + t1 * dir_vec_z;
                            // is inside the borders
                            if( intersection_z > z0_ - 0.5*z_ &&
                                intersection_z < z0_ + 0.5*z_    )
                            {
                                std::swap(distance.first,distance.second);
                                distance.first  =   t1;
                            }
                        }
                        else
                        {
                            intersection_z  =   particle->GetZ() + t2 * dir_vec_z;
                            // is inside the borders
                            if( intersection_z > z0_ - 0.5*z_ &&
                                intersection_z < z0_ + 0.5*z_    )
                            {
                                std::swap(distance.first,distance.second);
                                distance.first  =   t2;
                            }
                        }
                    }
                    // Now we have to check if the particle is on the border of
                    // the inner sphere
                    if(t1 == 0 )
                    {
                        // The particle is moving into the inner cylinder
                        if( t2 > 0 )
                        {
                            intersection_z  =   particle->GetZ() + t2 * dir_vec_z;
                            // is inside the borders
                            if( intersection_z > z0_ - 0.5*z_ &&
                                intersection_z < z0_ + 0.5*z_    )
                            {
                                std::swap(distance.first,distance.second);
                                distance.first  =   t2;
                            }

                        }
                        // if not we don't have to update distance.first
                    }
                    if(t2 == 0 )
                    {
                        // The particle is moving into the inner sphere
                        if( t1 > 0 )
                        {
                            intersection_z  =   particle->GetZ() + t1 * dir_vec_z;
                            // is inside the borders
                            if( intersection_z > z0_ - 0.5*z_ &&
                                intersection_z < z0_ + 0.5*z_    )
                            {
                                std::swap(distance.first,distance.second);
                                distance.first  =   t1;
                            }
                        }
                        // if not we don't have to update distance.first
                    }
                }
            }
        }
    }

    // Make a computer precision controll!
    // This is necessary cause due to numerical effects it meight be happen
    // that a particle which is located on a gemoetry border is treated as inside
    // or outside

    if( distance.first < GEOMETRY_PRECISION )
        distance.first  =   -1;
    if( distance.second < GEOMETRY_PRECISION )
        distance.second  =   -1;
    if( distance.first < 0 )
        std::swap(distance.first ,distance.second);

    return distance;
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

void Geometry::SetInnerRadius(double inner_radius) {
    inner_radius_ =   inner_radius;
}

void Geometry::SetRadius(double radius) {
    radius_ =   radius;
}

void Geometry::SetObject(string object) {
    object_ =   object;
}

void Geometry::SetHirarchy(unsigned int hirarchy) {
    hirarchy_ =   hirarchy;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//---------------------------------Destructor---------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Geometry::~Geometry(){}
