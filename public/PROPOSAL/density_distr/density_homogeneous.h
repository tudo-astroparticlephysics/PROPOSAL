#pragma once
#include "PROPOSAL/density_distr/density_distr.h"
 
namespace PROPOSAL{
 
class Density_homogeneous : public Density_distr
{
public:
    Density_homogeneous();
    Density_homogeneous(const Density_homogeneous&);
    Density_homogeneous(Vector3D fAxis, Vector3D fp0, std::function<double(double)> density_distribution);

    ~Density_homogeneous() {};

    Density_homogeneous* clone() const {return new Density_homogeneous(*this);};
 
    double Integrate(Vector3D xi, Vector3D direction, double res) const override;
    
};
 
}

