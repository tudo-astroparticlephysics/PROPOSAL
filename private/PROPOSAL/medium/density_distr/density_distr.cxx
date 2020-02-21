
#include <iostream>
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/medium/density_distr/density_distr.h"

using namespace PROPOSAL;

Density_distr::Density_distr() : axis_(CartesianAxis().clone()) {}

Density_distr::Density_distr(const Density_distr& density_distr)
    : axis_(density_distr.axis_->clone()) {}

Density_distr::Density_distr(const Axis& axis) : axis_(axis.clone()) {}

bool Density_distr::operator==(const Density_distr& dens_distr) const
{
    if (*axis_ != *dens_distr.axis_)
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

Axis::Axis(const Vector3D& fAxis, const Vector3D& fp0) : fAxis_(fAxis), fp0_(fp0) {}

Axis::Axis(const Axis& axis) : fAxis_(axis.fAxis_), fp0_(axis.fp0_) {}

bool Axis::operator==(const Axis& axis) const {
    if(fAxis_ != axis.fAxis_)
        return false;
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
    fp0_.SetCartesianCoordinates(0, 0, 0);
    fAxis_.SetSphericalCoordinates(1, 0, 0);
}

RadialAxis::RadialAxis(const Vector3D& fAxis, const Vector3D& fp0) : Axis(fAxis, fp0) {}

double RadialAxis::GetDepth(const Vector3D& xi) const {
    return (xi - fp0_).magnitude();
}

double RadialAxis::GetEffectiveDistance(const Vector3D& xi, const Vector3D& direction) const {
    Vector3D aux{xi - fp0_};
    aux.normalise();

    return -aux * direction;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%      Cartesian     %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CartesianAxis::CartesianAxis() : Axis() {
    fAxis_.SetCartesianCoordinates(1, 0, 0);
    fp0_.SetCartesianCoordinates(0, 0, 0);
}

CartesianAxis::CartesianAxis(const Vector3D& fAxis, const Vector3D& fp0) : Axis(fAxis, fp0) {}

double CartesianAxis::GetDepth(const Vector3D& xi) const {
    return fAxis_ * (xi - fp0_);
}

double CartesianAxis::GetEffectiveDistance(const Vector3D& xi,
                                           const Vector3D& direction) const {
    (void)xi;

    return fAxis_ * direction;
}
