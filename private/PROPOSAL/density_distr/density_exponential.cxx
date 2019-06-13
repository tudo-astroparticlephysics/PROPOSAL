#include "PROPOSAL/density_distr/density_exponential.h"
#include <cmath>
#include <iostream>


Density_exponential::Density_exponential(Vector3D fAxis, 
                                         Vector3D fp0, 
                                         double sigma):
    Density_distr(fAxis, fp0),
    sigma_(sigma)
{}

double Density_exponential::Correct(Vector3D xi, 
                                    Vector3D direction, 
                                    double res) const
{
    double phi = GetDepth(xi) / sigma_;
    double delta = GetAxis() * direction / sigma_;

    return 1. / delta * std::log( 1 + std::exp(-phi) * res * delta );
}

double Density_exponential::Integrate(Vector3D xi, 
                                      Vector3D direction, 
                                      double l) const
{
    double phi = GetDepth(xi) / sigma_;
    double delta = GetAxis() * direction / sigma_;

    return std::exp( phi + l * delta ) / delta;
}



double Density_exponential::Calculate(Vector3D xi, 
                                      Vector3D direction, 
                                      double distance) const
{
    return Integrate(xi, direction, distance) - Integrate(xi, direction, 0);
}
