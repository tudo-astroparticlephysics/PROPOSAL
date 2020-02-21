
#include <cmath>
#include <iostream>
#include "PROPOSAL/medium/density_distr/density_exponential.h"

using namespace PROPOSAL;

Density_exponential::Density_exponential(const Axis& axis, double sigma)
    : Density_distr(axis), sigma_(sigma) {}

double Density_exponential::GetDepth(const Vector3D& xi) const {
    return axis_->GetDepth(xi) / sigma_;
}

bool Density_exponential::compare(const Density_distr& dens_distr) const {
    const Density_exponential* dens_exp = dynamic_cast<const Density_exponential*>(&dens_distr);
    if(!dens_exp)
        return false;
    if(sigma_ != dens_exp->sigma_ )
        return false;
    return true;
}

double Density_exponential::GetEffectiveDistance(const Vector3D& xi,
                                                 const Vector3D& direction) const {
    return axis_->GetEffectiveDistance(xi, direction) / sigma_;
}

double Density_exponential::Correct(const Vector3D& xi,
                                    const Vector3D& direction,
                                    double res,
                                    double distance_to_border) const {
    (void)distance_to_border;

    double phi = GetDepth(xi);
    double delta = GetEffectiveDistance(xi, direction);

    double aux = 1. / delta * std::log(1 + std::exp(-phi) * res * delta);

    if (std::isnan(aux))
        throw DensityException("Next interaction point lies in infinite.");

    return aux;
}

double Density_exponential::Integrate(const Vector3D& xi,
                                      const Vector3D& direction,
                                      double l) const {
    double delta = GetEffectiveDistance(xi, direction);

    return std::exp(GetDepth(xi) + l * delta) / delta;
}

double Density_exponential::Calculate(const Vector3D& xi,
                                      const Vector3D& direction,
                                      double distance) const {
    return Integrate(xi, direction, distance) - Integrate(xi, direction, 0);
}

double Density_exponential::Evaluate(const Vector3D& xi) const {
    return std::exp(GetDepth(xi));
}
