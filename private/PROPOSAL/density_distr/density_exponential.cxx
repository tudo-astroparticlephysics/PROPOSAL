#include "PROPOSAL/density_distr/density_exponential.h"
#include <cmath>
#include <iostream>


Density_exponential::Density_exponential(const Axis& axis, double sigma):
    axis_(axis.clone()),
    sigma_(sigma)
{}

double Density_exponential::GetDepth(Vector3D xi) const 
{
    std::cout << "Depth: " << axis_->GetDepth(xi) / sigma_ << std::endl;
    return axis_->GetDepth(xi) / sigma_;
}

double Density_exponential::GetEffectiveDistance(Vector3D direction) const 
{
    std::cout << "Effective Distance: " << axis_->GetEffectiveDistance(direction) / sigma_ << std::endl;
    return axis_->GetEffectiveDistance(direction) / sigma_;
}

double Density_exponential::Correct(Vector3D xi, 
                                    Vector3D direction, 
                                    double res) const 
{
    double phi = GetDepth(xi);
    double delta = GetEffectiveDistance(direction);

    return 1. / delta * std::log( 1 + std::exp(-phi) * res * delta );
}

double Density_exponential::Integrate(Vector3D xi, 
                                      Vector3D direction, 
                                      double l) const 
{
    double phi = GetDepth(xi);
    double delta = GetEffectiveDistance(direction);

    return std::exp( phi + l * delta ) / delta;
}



double Density_exponential::Calculate(Vector3D xi, 
                                      Vector3D direction, 
                                      double distance) const 
{
    return Integrate(xi, direction, distance) - Integrate(xi, direction, 0);
}
