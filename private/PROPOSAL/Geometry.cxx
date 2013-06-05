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
    double distance =   0;
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
