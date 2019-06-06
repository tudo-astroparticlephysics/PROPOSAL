#pragma once
#include "PROPOSAL/density_distr/density_distr.h"
 
namespace PROPOSAL{
 
class Density_homogeneous : public Density_distr
{
public:
    Density_homogeneous(double rho);
 
    double Integrate(double x_i, double res);
};
 
}

