#include "PROPOSAL/density_distr/density_distr.h"
#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/Logging.h"
#include <iostream>

Density_distr::Density_distr():
    axis_(CartesianAxis().clone())
{}

Density_distr::Density_distr(const Density_distr& density_distr):
    axis_(density_distr.axis_)
{}
                                                    
Density_distr::Density_distr(const Axis& axis):
    axis_(axis.clone())
{}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%        Axis        %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Axis::Axis()
{}

Axis::Axis(Vector3D fAxis, Vector3D fp0):
    fAxis_(fAxis),
    fp0_(fp0)
{}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%       Radial       %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RadialAxis::RadialAxis():
    Axis() 
{
    fAxis_.SetSphericalCoordinates(1,0,0);
}

RadialAxis::RadialAxis(Vector3D fp0):
    Axis()
{
    fp0_ = fp0;
    fAxis_.SetSphericalCoordinates(1,0,0);
}

double RadialAxis::GetDepth(Vector3D xi) const
{
    // angle dependencies require more precise investigations MaxSac (2019-06-18)
    xi.CalculateSphericalCoordinates();

    return xi.GetRadius() - fp0_.GetRadius();
}

double RadialAxis::GetEffectiveDistance(Vector3D direction) const 
{
    direction.CalculateSphericalCoordinates();

    return (fAxis_ - direction) * direction; 

    return fAxis_.GetRadius() * direction.GetRadius();
}


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%      Cartesian     %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CartesianAxis::CartesianAxis():
    Axis()
{
    fAxis_.SetCartesianCoordinates(1,0,0); 
    fp0_.SetCartesianCoordinates(0,0,0);
}

CartesianAxis::CartesianAxis(Vector3D fAxis, Vector3D fp0):
    Axis(fAxis,fp0)
{}

double CartesianAxis::GetDepth(Vector3D xi) const 
{
    return fAxis_ * (xi - fp0_);
}

double CartesianAxis::GetEffectiveDistance(Vector3D direction) const 
{
    return fAxis_ * direction;
}
