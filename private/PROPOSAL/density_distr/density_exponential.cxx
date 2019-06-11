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

    return 1. / delta * ( std::log( std::exp(phi) + res * delta ) - phi );
}

double Density_exponential::Integrate(Vector3D xr, 
                                      Vector3D direction, 
                                      double l) const
{
    double aux1 = ( direction * fAxis_ ) / sigma_;
    double aux2 = GetDepth(xr) / sigma_;

    return std::exp( aux1 * l + aux2 ) / aux1;
}



double Density_exponential::Calculate(Vector3D xi, 
                                      Vector3D direction, 
                                      double distance) const
{
    /* Pruefe ob es sinn macht das hier ein fabs steht. */
    return std::fabs(Integrate(xi, direction, 0) - Integrate(xi, direction, distance));
}
