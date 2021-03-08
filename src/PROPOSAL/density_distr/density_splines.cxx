
#include <functional>
#include "PROPOSAL/density_distr/density_splines.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/math/Cartesian3D.h"
using namespace PROPOSAL;

Density_splines::Density_splines(const Axis& axis, const Spline& splines, double massDensity)
    : Density_distr(axis, massDensity),
      spline_(splines.clone()),
      integrated_spline_(splines.clone()) {
    integrated_spline_->Antiderivative(0);
}

Density_splines::Density_splines(const PROPOSAL::Axis& axis, const PROPOSAL::Spline& splines, const Medium& medium)
    : Density_splines(axis, splines, medium.GetMassDensity()) {}

Density_splines::Density_splines(const Density_splines& dens_splines)
    : Density_distr(dens_splines),
      spline_(dens_splines.spline_),
      integrated_spline_(dens_splines.integrated_spline_) {}

Density_splines::Density_splines(const nlohmann::json& config) : Density_distr(config) {
    if(!config.contains("spline_type"))
        throw std::invalid_argument("Density_splines: Type of splines must be specified using spline_type");
    std::string spline_type = config["spline_type"];
    if(spline_type == "linear") {
        spline_ = new Linear_Spline(config);
        integrated_spline_ = new Linear_Spline(config);
    } else if (spline_type == "cubic") {
        spline_ = new Cubic_Spline(config);
        integrated_spline_ = new Cubic_Spline(config);
    } else {
        throw std::invalid_argument("Density_splines: Type of spline must be linear or cubic");
    }
    integrated_spline_->Antiderivative(0);
}


bool Density_splines::compare(const Density_distr& dens_distr) const {
    const Density_splines* dens_splines= dynamic_cast<const Density_splines*>(&dens_distr);
    if(!dens_splines)
        return false;
    if( spline_!= dens_splines->spline_)
        return false;
    if( integrated_spline_!= dens_splines->integrated_spline_)
        return false;
    return true;
}

double Density_splines::Helper_function(const Vector3D& xi,
                                        const Vector3D& direction,
                                        double res,
                                        double l) const {
    return Integrate(xi, direction, 0) - Integrate(xi, direction, l) + res;
}

double Density_splines::helper_function(const Vector3D& xi,
                                        const Vector3D& direction,
                                        double res,
                                        double l) const {
    (void)res;
    auto dir_cartesian = Cartesian3D(direction);
    return Evaluate(xi) - Evaluate(xi + l * dir_cartesian);
}

double Density_splines::Correct(const Vector3D& xi,
                                const Vector3D& direction,
                                double res,
                                double distance_to_border) const {
    auto F = [&](double l) {return Helper_function(xi, direction, res, l);};
    auto dF = [&](double l) {return helper_function(xi, direction, res, l);};

    try {
        res =
            NewtonRaphson(F, dF, 0, distance_to_border, distance_to_border / 2);
    } catch (MathException& e) {
        throw DensityException("Next interaction point lies in infinite.");
    }

    return res;
}

double Density_splines::Integrate(const Vector3D& xi,
                                  const Vector3D& direction,
                                  double l) const {
    double delta = axis_->GetEffectiveDistance(xi, direction);

    return massDensity_ / (delta * delta) *
           integrated_spline_->evaluate(axis_->GetDepth(xi) + l * delta);
}

double Density_splines::Calculate(const Vector3D& xi,
                                  const Vector3D& direction,
                                  double distance) const {
    return Integrate(xi, direction, distance) - Integrate(xi, direction, 0);
}

double Density_splines::Evaluate(const Vector3D& xi) const {
    return massDensity_ * spline_->evaluate(axis_->GetDepth(xi));
}
