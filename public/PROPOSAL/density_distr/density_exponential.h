#pragma once
 
#include "PROPOSAL/density_distr/density_distr.h"
#include "PROPOSAL/math/Vector3D.h"
 
 
class Density_exponential: public Density_distr
{
public:
    Density_exponential(Vector3D fAxis, Vector3D fp0, std::function<double(double)> density_distribution, double sigma);
 
    double Integrate(Vector3D xi, Vector3D direction, double res);

private:
    double sigma_;
};


