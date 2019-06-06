#include "PROPOSAL/density_distr/density_homogeneous.h"                                                                                                    
 
Density_homogeneous::Density_homogeneous(double rho):
    Dens_distr(rho)
{};
 
double Density_homogeneous::Integrate(double x_i, double res):
    x_i_(x_i),
    res_(res)
{
    return x_i_ - res_ / rho_ 
}

