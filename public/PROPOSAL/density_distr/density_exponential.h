#pragma once
 
#include "PROPOSAL/density_distr/density_distr.h"
#include "PROPOSAL/math/Vector3D.h"
 
 
class Density_exponential: public Density_distr
{
public:
    Density_exponential(Vector3D fAxis, Vector3D fp0, std::function<double(double)> density_distribution, double sigma);
 
    ~Density_exponential() {};

    Density_exponential* clone() const {return new Density_exponential(*this);};

    double Integrate(Vector3D xi, Vector3D direction, double res) const;

private:
    double sigma_;
};


