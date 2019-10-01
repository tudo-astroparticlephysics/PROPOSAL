
#include "PROPOSAL/medium/density_distr/density_homogeneous.h"
#include <cmath>
#include <iostream>

using namespace PROPOSAL;

Density_homogeneous::Density_homogeneous()
    : Density_distr(), correction_factor_(1.0) {}

Density_homogeneous::Density_homogeneous(double correction_factor)
    : Density_distr(), correction_factor_(correction_factor) {}

Density_homogeneous::Density_homogeneous(const Density_homogeneous& dens_distr)
    : Density_distr(dens_distr),
      correction_factor_(dens_distr.correction_factor_) {}

double Density_homogeneous::Correct(Vector3D xi,
                                    Vector3D direction,
                                    double res,
                                    double distance_to_border) const {
    (void)xi;
    (void)direction;
    (void)distance_to_border;

    return res / correction_factor_;
}

double Density_homogeneous::Integrate(Vector3D xi,
                                      Vector3D direction,
                                      double l) const {
    (void)xi;
    (void)direction;

    return correction_factor_ * l;
}

double Density_homogeneous::Calculate(Vector3D xi,
                                      Vector3D direction,
                                      double distance) const {
    return Integrate(xi, direction, distance) - Integrate(xi, direction, 0);
}

double Density_homogeneous::Evaluate(Vector3D xi) const {
    (void)xi;

    return correction_factor_;
}
