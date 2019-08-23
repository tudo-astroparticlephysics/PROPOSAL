#include "PROPOSAL/density_distr/density_splines.h"
#include <algorithm>
#include <functional>
#include <iostream>
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/math/Spline.h"

Density_splines::Density_splines(const Axis& axis, const Spline& splines)
    : Density_distr(axis),
      spline_(splines.clone()),
      integrated_spline_(splines.clone()) {
    integrated_spline_->Antiderivative(0);
}

Density_splines::Density_splines(const Density_splines& dens_splines)
    : Density_distr(dens_splines),
      spline_(dens_splines.spline_),
      integrated_spline_(dens_splines.integrated_spline_) {}

double Density_splines::Helper_function(Vector3D xi,
                                        Vector3D direction,
                                        double res,
                                        double l) const {
    return Integrate(xi, direction, 0) - Integrate(xi, direction, l) + res;
}

double Density_splines::helper_function(Vector3D xi,
                                        Vector3D direction,
                                        double res,
                                        double l) const {
    return Evaluate(xi, direction, 0) - Evaluate(xi, direction, l);
}

double Density_splines::Correct(Vector3D xi,
                                Vector3D direction,
                                double res,
                                double distance_to_border) const {
    std::function<double(double)> F =
        std::bind(&Density_splines::Helper_function, this, xi, direction, res,
                  std::placeholders::_1);

    std::function<double(double)> dF =
        std::bind(&Density_splines::helper_function, this, xi, direction, res,
                  std::placeholders::_1);

    try {
        res =
            NewtonRaphson(F, dF, 0, distance_to_border, distance_to_border / 2);
    } catch (MathException& e) {
        throw DensityException("Next interaction point lies in infinite.");
    }

    return res;
}

double Density_splines::Integrate(Vector3D xi,
                                  Vector3D direction,
                                  double l) const {
    double delta = axis_->GetEffectiveDistance(xi, direction);

    return 1 / (delta * delta) *
           integrated_spline_->evaluate(axis_->GetDepth(xi) + l * delta);
}

double Density_splines::Evaluate(Vector3D xi,
                                 Vector3D direction,
                                 double l) const {
    double delta = axis_->GetEffectiveDistance(xi, direction);

    return spline_->evaluate(axis_->GetDepth(xi) + l * delta);
}

double Density_splines::Calculate(Vector3D xi,
                                  Vector3D direction,
                                  double distance) const {
    // return Integrate(xi, direction, distance) - Integrate(xi, direction, 0);
    double aux =
        Integrate(xi, direction, distance) - Integrate(xi, direction, 0);
    std::cout << "Calculate(" << distance << "): " << aux << std::endl;
    return aux;
}

