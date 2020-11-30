
#include <cmath>
#include <iostream>
#include "PROPOSAL/density_distr/density_homogeneous.h"
#include "PROPOSAL/medium/Medium.h"
using namespace PROPOSAL;

Density_homogeneous::Density_homogeneous()
    : Density_distr(1.), correction_factor_(1.0) {}

Density_homogeneous::Density_homogeneous(double massDensity, double correction_factor)
    : Density_distr(massDensity), correction_factor_(correction_factor) {}

Density_homogeneous::Density_homogeneous(const Medium &medium, double correction_factor)
    : Density_homogeneous(medium.GetMassDensity(), correction_factor) {}

Density_homogeneous::Density_homogeneous(const Density_homogeneous& dens_distr)
    : Density_distr(dens_distr),
      correction_factor_(dens_distr.correction_factor_) {}

Density_homogeneous::Density_homogeneous(const nlohmann::json& config) : Density_distr() {
    correction_factor_ = config.value("correction_factor", 1.);
    if (config.contains("massDensity")) {
        assert(config["massDensity"].is_number());
        massDensity_ = config["massDensity"].get<double>();
    } else {
        throw std::invalid_argument("Density_distr: MassDensity must be defined in json");
    }
}


bool Density_homogeneous::compare(const Density_distr& dens_distr) const {
    const Density_homogeneous* dens_homogen = dynamic_cast<const Density_homogeneous*>(&dens_distr);
    if(!dens_homogen)
        return false;
    if(correction_factor_ != dens_homogen->correction_factor_ )
        return false;
    return true;
}

double Density_homogeneous::Correct(const Vector3D& xi,
                                    const Vector3D& direction,
                                    double res,
                                    double distance_to_border) const {
    (void)xi;
    (void)direction;
    (void)distance_to_border;

    return res / (correction_factor_ * massDensity_);
}

double Density_homogeneous::Integrate(const Vector3D& xi,
                                      const Vector3D& direction,
                                      double l) const {
    (void)xi;
    (void)direction;

    return correction_factor_ * massDensity_ * l;
}

double Density_homogeneous::Calculate(const Vector3D& xi,
                                      const Vector3D& direction,
                                      double distance) const {
    return Integrate(xi, direction, distance) - Integrate(xi, direction, 0);
}

double Density_homogeneous::Evaluate(const Vector3D& xi) const {
    (void)xi;

    return correction_factor_ * massDensity_;
}
