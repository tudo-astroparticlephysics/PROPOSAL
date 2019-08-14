#include "PROPOSAL/density_distr/density_exponential.h"
#include <cmath>
#include <iostream>


Density_exponential::Density_exponential(const Axis& axis, double sigma):
    Density_distr(axis),
    sigma_(sigma)
{}

double Density_exponential::GetDepth(Vector3D xi) const 
{
    return axis_->GetDepth(xi) / sigma_;
}

double Density_exponential::GetEffectiveDistance(Vector3D xi, Vector3D direction) const 
{
    return axis_->GetEffectiveDistance(xi, direction) / sigma_;
}

double Density_exponential::Correct(Vector3D xi, 
                                    Vector3D direction, 
                                    double res,
                                    double distance_to_border) const 
{
    double phi = GetDepth(xi);
    double delta = GetEffectiveDistance(xi, direction);

    return 1. / delta * std::log( 1 + std::exp(-phi) * res * delta );
}

double Density_exponential::Integrate(Vector3D xi, 
                                      Vector3D direction, 
                                      double l) const 
{
    double delta = GetEffectiveDistance(xi, direction);

    return std::exp( GetDepth(xi) + l * delta ) / delta;
}



double Density_exponential::Calculate(Vector3D xi, 
                                      Vector3D direction, 
                                      double distance) const 
{
    return Integrate(xi, direction, distance) - Integrate(xi, direction, 0);
}

double Density_exponential::GetCorrection(Vector3D xi) const
{
    return std::exp(GetDepth(xi));
}
