#pragma once

#include "PROPOSAL/density_distr/density_distr.h"

namespace PROPOSAL{

class Density_exponential: public Density_distr
{
public:
    Density_exponential(double rho, double sigma, double mu);

    double Integrate(double x_i, double res);

private:
    double sigma_;
    double mu_;
}

}


