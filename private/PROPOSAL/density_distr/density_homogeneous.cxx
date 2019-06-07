#include "PROPOSAL/density_distr/density_homogeneous.h"                                                                                                    
 
Density_homogeneous::Density_homogeneous():
    Density_distr()
{}

Density_homogeneous::Density_homogeneous(const Density_homogeneous& density_homogeneous):
    Density_distr(density_homogeneous.fAxis_, density_homogeneous.fp0_, density_homogeneous.density_distribution_)
{}

Density_homogeneous::Density_homogeneous(Vector3D fAxis, 
                                         Vector3D fp0, 
                                         std::function<double(double)> density_distribution):
    Density_distr(fAxis, fp0, density_distribution)
{}
 
double Density_homogeneous::Integrate(Vector3D xi, Vector3D direction, double res) const
{
    (void) direction;

    return xi.magnitude() - res; 
}

