#include "PROPOSAL/density_distr/density_homogeneous.h"                                                                                                    
 
Density_homogeneous::Density_homogeneous():
    Density_distr()
{ }

Density_homogeneous::Density_homogeneous(Vector3D fAxis, 
                                         Vector3D fp0, 
                                         std::function<double(double)> density_distribution):
    Density_distr(fAxis, fp0, density_distribution)
{}
 
double Density_homogeneous::Integrate(Vector3D xi, Vector3D direction, double res)
{
    (void) direction;

    return xi.magnitude() - res; 
}

