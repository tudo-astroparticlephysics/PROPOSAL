
#include <cmath>
#include <iostream>
#include "PROPOSAL/density_distr/density_homogeneous.h"
#include "PROPOSAL/medium/Medium.h"
using namespace PROPOSAL;

Density_homogeneous::Density_homogeneous(double massDensity)
    : Density_distr(massDensity) {}

Density_homogeneous::Density_homogeneous(const Medium &medium, double correction_factor)
    : Density_homogeneous(medium.GetMassDensity() * correction_factor) {}

Density_homogeneous::Density_homogeneous(const Density_homogeneous& dens_distr)
    : Density_distr(dens_distr) {}

Density_homogeneous::Density_homogeneous(const nlohmann::json& config) : Density_distr() {
    if (config.contains("mass_density")) {
        assert(config["mass_density"].is_number());
        massDensity_ = config["mass_density"].get<double>();
    } else {
        throw std::invalid_argument("Density_distr: mass_density must be defined in json");
    }
}

bool Density_homogeneous::compare(const Density_distr& dens_distr) const {
    const Density_homogeneous* dens_homogen = dynamic_cast<const Density_homogeneous*>(&dens_distr);
    if(!dens_homogen)
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

    return res / massDensity_;
}

double Density_homogeneous::Integrate(const Vector3D& xi,
                                      const Vector3D& direction,
                                      double l) const {
    (void)xi;
    (void)direction;

    return massDensity_ * l;
}

double Density_homogeneous::Calculate(const Vector3D& xi,
                                      const Vector3D& direction,
                                      double distance) const {
    return Integrate(xi, direction, distance) - Integrate(xi, direction, 0);
}

double Density_homogeneous::Evaluate(const Vector3D& xi) const {
    (void)xi;

    return massDensity_;
}
