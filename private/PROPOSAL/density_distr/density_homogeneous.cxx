#include "PROPOSAL/density_distr/density_homogeneous.h"                                                                                                    
#include <cmath>
#include <iostream>
 
Density_homogeneous::Density_homogeneous():
    Density_distr()
{}

Density_homogeneous::Density_homogeneous(const Density_homogeneous& density_homogeneous):
    Density_distr(density_homogeneous.fAxis_, density_homogeneous.fp0_)
{}

Density_homogeneous::Density_homogeneous(Vector3D fAxis, Vector3D fp0):
    Density_distr(fAxis, fp0)
{}
 
double Density_homogeneous::Correct(Vector3D xi, 
                                    Vector3D direction, 
                                    double res) const
{
    (void) xi;
    (void) direction;

    return res; 
}

double Density_homogeneous::Integrate(Vector3D xi, 
                                      Vector3D direction, 
                                      double distance) const
{
    (void) xi;
    (void) direction;

    return distance; 
}

double Density_homogeneous::Calculate(Vector3D xi, 
                                      Vector3D direction, 
                                      double distance) const 
{
    (void)xi;
    (void)direction;

    /* std::cout << distance << std::endl; */
    return distance;
}
