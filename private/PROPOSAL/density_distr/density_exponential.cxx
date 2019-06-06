#include "PROPOSAL/density_distr/density_exponential.h"
#include <cmath>


Density_exponential::Density_exponential(Vector3D fAxis, 
                                         Vector3D fp0, 
                                         std::function<double(double)> density_distribution, 
                                         double sigma):
    Density_distr(fAxis, fp0, density_distribution),
    sigma_(sigma)
{
}

double Density_exponential::Integrate(Vector3D xi, Vector3D direction, double res)
{
    double aux1 = GetDepth(xi) / sigma_;
    double aux2 = fAxis_ * direction / sigma_;

    return 1 / aux2 * std::log( std::exp(aux1) - res * aux2 - aux1 );
}



