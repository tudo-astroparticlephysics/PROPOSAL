#pragma once
#include "PROPOSAL/density_distr/density_distr.h"
 
namespace PROPOSAL{
 
class Density_homogeneous : public Density_distr
{
public:
    Density_homogeneous();
    Density_homogeneous(const Density_homogeneous&);
    Density_homogeneous(const Axis& axis);

    ~Density_homogeneous() {};

    Density_homogeneous* clone() const {return new Density_homogeneous(*this);};
 
    double Correct(Vector3D xi, Vector3D direction, double res, double distance_to_border) const override;
    double Integrate(Vector3D xi, Vector3D direction, double res) const override;
    double Calculate(Vector3D xi, Vector3D direction, double distance) const override;
    double GetCorrection(Vector3D xi) const override;
};
 
}

