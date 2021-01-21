
#include <functional>
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/density_distr/density_polynomial.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/math/Cartesian3D.h"
using namespace PROPOSAL;
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%% Polynomial-Density %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Density_polynomial::Density_polynomial(const Axis& axis, const Polynom& polynom, double massDensity)
    : Density_distr(axis, massDensity),
      polynom_(polynom),
      Polynom_(polynom_.GetAntiderivative(0)),
      density_distribution(polynom_.GetFunction()),
      antiderived_density_distribution(Polynom_.GetFunction()) {}

Density_polynomial::Density_polynomial(const Axis& axis, const Polynom& polynom, const Medium& medium)
    : Density_polynomial(axis, polynom, medium.GetMassDensity()) {}

Density_polynomial::Density_polynomial(const Density_polynomial& dens)
    : Density_distr(dens),
      polynom_(dens.polynom_),
      Polynom_(dens.Polynom_),
      density_distribution(polynom_.GetFunction()),
      antiderived_density_distribution(Polynom_.GetFunction()) {}

Density_polynomial::Density_polynomial(const nlohmann::json& config) : Density_distr(config),
    polynom_(config),
    Polynom_(polynom_.GetAntiderivative(0)),
    density_distribution(polynom_.GetFunction()),
    antiderived_density_distribution(Polynom_.GetFunction()) {}


Density_polynomial::~Density_polynomial() {}

bool Density_polynomial::compare(const Density_distr& dens_distr) const {
    const Density_polynomial* dens_poly = dynamic_cast<const Density_polynomial*>(&dens_distr);
    if(!dens_poly)
        return false;
    if( polynom_ != dens_poly->polynom_ )
        return false;
    return true;
}

double Density_polynomial::Helper_function(const Vector3D& xi,
                                           const Vector3D& direction,
                                           double res,
                                           double l) const {
    return Integrate(xi, direction, 0) - Integrate(xi, direction, l) + res;
}

double Density_polynomial::helper_function(const Vector3D& xi,
                                           const Vector3D& direction,
                                           double res,
                                           double l) const {
    (void)res;
    Cartesian3D dir_cartesian(direction);
    return Evaluate(xi) - Evaluate(xi + l * dir_cartesian);
}

double Density_polynomial::Correct(const Vector3D& xi,
                                   const Vector3D& direction,
                                   double res,
                                   double distance_to_border) const {
    auto F = [&](double l) {return Helper_function(xi, direction, res, l);};
    auto dF = [&](double l) {return helper_function(xi, direction, res, l);};

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

double Density_polynomial::Integrate(const Vector3D& xi,
                                     const Vector3D& direction,
                                     double l) const {
    double delta = axis_->GetEffectiveDistance(xi, direction);

    return massDensity_ * antiderived_density_distribution(axis_->GetDepth(xi) + l * delta) /
           (delta * delta);
}

double Density_polynomial::Calculate(const Vector3D& xi,
                                     const Vector3D& direction,
                                     double distance) const {
    return Integrate(xi, direction, distance) - Integrate(xi, direction, 0);
}

double Density_polynomial::Evaluate(const Vector3D& xi) const {
    return massDensity_ * density_distribution(axis_->GetDepth(xi));
}
