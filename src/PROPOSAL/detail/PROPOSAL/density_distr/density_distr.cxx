
#include <iostream>
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/density_distr/density_distr.h"
#include "PROPOSAL/density_distr/density_exponential.h"
#include "PROPOSAL/density_distr/density_homogeneous.h"
#include "PROPOSAL/density_distr/density_polynomial.h"
#include "PROPOSAL/density_distr/density_splines.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/math/Cartesian3D.h"

using namespace PROPOSAL;

Density_distr::Density_distr() : axis_(CartesianAxis().clone()), massDensity_(1.) {}

Density_distr::Density_distr(const Density_distr& density_distr)
    : axis_(density_distr.axis_->clone()), massDensity_(density_distr.massDensity_) {}

Density_distr::Density_distr(const nlohmann::json& config) {
    if(!config.contains("axis_type"))
        throw std::invalid_argument("Axis: Type of axis must be specified using axis_type");
    std::string axis_type = config["axis_type"];
    if(axis_type == "radial") {
        axis_ = new RadialAxis(config);
    } else if (axis_type == "cartesian") {
        axis_ = new CartesianAxis(config);
    } else {
        throw std::invalid_argument("Axis: Type of axis must be radial or cartesian");
    }
    if (config.contains("mass_density")) {
        assert(config["mass_density"].is_number());
        massDensity_ = config["mass_density"].get<double>();
    } else {
        throw std::invalid_argument("Density_distr: mass_density must be defined in json");
    }
}

Density_distr::Density_distr(double massDensity) : axis_(CartesianAxis().clone()), massDensity_(massDensity) {}

Density_distr::Density_distr(const Axis& axis, double massDensity) : axis_(axis.clone()), massDensity_(massDensity) {}

bool Density_distr::operator==(const Density_distr& dens_distr) const
{
    if (*axis_ != *dens_distr.axis_)
        return false;
    if (massDensity_ != dens_distr.massDensity_)
        return false;
    if (!this->compare(dens_distr) )
        return false;
    return true;
}


bool Density_distr::operator!=(const Density_distr& dens_distr) const {
    return !(*this == dens_distr);
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%        Axis        %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Axis::Axis() {}

Axis::Axis(const Vector3D& fp0) : fp0_(fp0) {}

Axis::Axis(const nlohmann::json& config) {
    if(!config.contains("fp0")) throw std::invalid_argument("Axis: No fp0 found.");
    fp0_ = Cartesian3D(config.at("fp0"));
}


Axis::Axis(const Axis& axis) : fp0_(axis.fp0_) {}

bool Axis::operator==(const Axis& axis) const {
    if(fp0_ != axis.fp0_)
        return false;
    return true;
}

bool Axis::operator!=(const Axis& axis) const {
    return !(*this == axis);
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%       Radial       %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RadialAxis::RadialAxis() : Axis() {
    fp0_.SetCoordinates({0, 0, 0});
}

RadialAxis::RadialAxis(const Vector3D& fp0) : Axis(fp0) {}

RadialAxis::RadialAxis(const nlohmann::json& config) : Axis(config) {}

double RadialAxis::GetDepth(const Vector3D& xi) const {
    return (xi - fp0_).magnitude();
}

double RadialAxis::GetEffectiveDistance(const Vector3D& xi, const Vector3D& direction) const {
    Cartesian3D aux{Cartesian3D(xi) - fp0_};
    aux.normalize();

    return -aux * direction;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%      Cartesian     %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CartesianAxis::CartesianAxis() : Axis() {
    fAxis_.SetCoordinates({1, 0, 0});
    fp0_.SetCoordinates({0, 0, 0});
}

CartesianAxis::CartesianAxis(const Vector3D& fAxis, const Vector3D& fp0) :  Axis(fp0) , fAxis_(fAxis) {}

CartesianAxis::CartesianAxis(const nlohmann::json& config) : Axis(config) {
    if(!config.contains("fAxis")) throw std::invalid_argument("Axis: No fAxis found.");
    fAxis_ = Cartesian3D(config.at("fAxis"));
}

bool CartesianAxis::operator==(const CartesianAxis& axis) const {
    if(fAxis_ != axis.fAxis_)
        return false;
    return Axis::operator==(axis);
}

bool CartesianAxis::operator!=(const CartesianAxis& axis) const {
    return !(*this == axis);
}

double CartesianAxis::GetDepth(const Vector3D& xi) const {
    return fAxis_ * (xi - fp0_);
}

double CartesianAxis::GetEffectiveDistance(const Vector3D& xi,
                                           const Vector3D& direction) const {
    (void)xi;

    return fAxis_ * direction;
}

std::shared_ptr<Density_distr> PROPOSAL::CreateDensityDistribution(const nlohmann::json& config) {
    if (config.contains("type")){
        std::string density_distr_type = config["type"];
        if (density_distr_type == "exponential") {
            return std::make_shared<Density_exponential>(config);
        } else if (density_distr_type == "homogeneous") {
            return std::make_shared<Density_homogeneous>(config);
        } else if (density_distr_type == "polynomial") {
            return std::make_shared<Density_polynomial>(config);
        } else if (density_distr_type == "spline") {
            return std::make_shared<Density_splines>(config);
        } else {
            throw std::invalid_argument("Density distribution config file must contain a paremeter called "
                                        "'type' with one of the keywords 'exponential', 'homogeneous',"
                                        " 'polynomial' or 'spline'.");
        }
    }
    else{
        throw std::invalid_argument("Density distribution config file must contain a parameter called 'type'.");
    }
}
