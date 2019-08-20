#include "PROPOSAL/density_distr/density_homogeneous.h"
#include <cmath>
#include <iostream>

Density_homogeneous::Density_homogeneous() : Density_distr() {}

Density_homogeneous::Density_homogeneous(const Density_homogeneous& dens_distr)
    : Density_distr() {
  (void)dens_distr;
}

double Density_homogeneous::Correct(Vector3D xi, Vector3D direction, double res,
                                    double distance_to_border) const {
  (void)xi;
  (void)direction;

  return res;
}

double Density_homogeneous::Integrate(Vector3D xi, Vector3D direction,
                                      double distance) const {
  (void)xi;
  (void)direction;

  return distance;
}

double Density_homogeneous::Calculate(Vector3D xi, Vector3D direction,
                                      double distance) const {
  (void)xi;
  (void)direction;

  return distance;
}

double Density_homogeneous::Evaluate(Vector3D xi, Vector3D direction,
                                     double l) const {
  (void)xi;
  (void)direction;
  (void)l;

  return 1;
}
