#include "PROPOSAL/density_distr/density_polynomial.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include "PROPOSAL/math/Function.h"
#include "PROPOSAL/math/MathMethods.h"

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%% Polynomial-Density %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Density_polynomial::Density_polynomial(const Axis& axis, const Polynom& polynom)
    : Density_distr(axis),
      polynom_(polynom),
      Polynom_(polynom_.GetAntiderivative(0)),
      density_distribution(polynom_.GetFunction()),
      antiderived_density_distribution(Polynom_.GetFunction()) {}

Density_polynomial::Density_polynomial(const Density_polynomial& dens)
    : Density_distr(dens),
      polynom_(dens.polynom_),
      Polynom_(dens.Polynom_),
      density_distribution(polynom_.GetFunction()),
      antiderived_density_distribution(Polynom_.GetFunction()) {}

double Density_polynomial::Helper_function(Vector3D xi,
                                           Vector3D direction,
                                           double res,
                                           double l) const {
    return Integrate(xi, direction, 0) - Integrate(xi, direction, l) + res;
}

double Density_polynomial::helper_function(Vector3D xi,
                                           Vector3D direction,
                                           double res,
                                           double l) const {
    return Evaluate(xi, direction, 0) - Evaluate(xi, direction, l);
}

double Density_polynomial::Correct(Vector3D xi,
                                   Vector3D direction,
                                   double res,
                                   double distance_to_border) const {
    std::function<double(double)> F =
        std::bind(&Density_polynomial::Helper_function, this, xi, direction,
                  res, std::placeholders::_1);

    std::function<double(double)> dF =
        std::bind(&Density_polynomial::helper_function, this, xi, direction,
                  res, std::placeholders::_1);

    // check if direction * axis larger or less than zero
    // direction * fAxis_

    try {
        res =
            NewtonRaphson(F, dF, 0, distance_to_border, distance_to_border / 2);
    } catch (MathException& e) {
        throw DensityException("Next interaction point lies in infinite.");
    }

    return res;
}

double Density_polynomial::Integrate(Vector3D xi,
                                     Vector3D direction,
                                     double l) const {
    double delta = axis_->GetEffectiveDistance(xi, direction);

    return 1 / (delta * delta) *
           antiderived_density_distribution(axis_->GetDepth(xi) + l * delta);
}

double Density_polynomial::Evaluate(Vector3D xi,
                                    Vector3D direction,
                                    double l) const {
    double delta = axis_->GetEffectiveDistance(xi, direction);

    return density_distribution(axis_->GetDepth(xi) + l * delta);
}

double Density_polynomial::Calculate(Vector3D xi,
                                     Vector3D direction,
                                     double distance) const {
    // return Integrate(xi, direction, distance) - Integrate(xi, direction, 0);
    double aux =
        Integrate(xi, direction, distance) - Integrate(xi, direction, 0);
    std::cout << "Calculate(" << distance << "): " << aux << std::endl;
    return aux;
}
