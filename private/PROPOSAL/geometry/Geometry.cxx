/*
 * Geometry.cxx
 *
 *  Created on: 05.06.2013
 *      Author: koehne
 */

#include <cmath>
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/geometry/Geometry.h"

#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

/******************************************************************************
 *                                  OStream                                    *
 ******************************************************************************/

namespace PROPOSAL {

std::ostream& operator<<(std::ostream& os, Geometry const& geometry)
{
    std::stringstream ss;
    ss << " Geometry (" << &geometry << ") ";
    os << Helper::Centered(60, ss.str()) << '\n';

    os << geometry.name_ << std::endl;
    os << "Position:\n" << geometry.position_ << '\n';
    os << "Hierarchy:\t" << geometry.hierarchy_ << '\n';

    geometry.print(os);

    os << Helper::Centered(60, "");

    return os;
}

} // namespace PROPOSAL

/******************************************************************************
 *                                  Geometry                                   *
 ******************************************************************************/

Geometry::Geometry(const std::string name)
    : position_(Vector3D())
    , name_(name)
    , hierarchy_(0)
{
}

Geometry::Geometry(const std::string name, const Vector3D position)
    : position_(100. * position)
    , name_(name)
    , hierarchy_(0)
{
}

Geometry::Geometry(const Geometry& geometry)
    : position_(geometry.position_)
    , name_(geometry.name_)
    , hierarchy_(geometry.hierarchy_)
{
}

// ------------------------------------------------------------------------- //
void Geometry::swap(Geometry& geometry)
{
    position_.swap(geometry.position_);
    name_.swap(geometry.name_);
    std::swap(hierarchy_, geometry.hierarchy_);
}

// ------------------------------------------------------------------------- //
Geometry& Geometry::operator=(const Geometry& geometry)
{
    if (this != &geometry)
    {
        position_ = geometry.position_;
        name_     = geometry.name_;
        hierarchy_ = geometry.hierarchy_;
    }

    return *this;
}

// ------------------------------------------------------------------------- //
bool Geometry::operator==(const Geometry& geometry) const
{
    if (position_ != geometry.position_)
        return false;
    else if (name_.compare(geometry.name_) != 0)
        return false;
    else if (hierarchy_ != geometry.hierarchy_)
        return false;
    else
        return this->compare(geometry);
}

// ------------------------------------------------------------------------- //
bool Geometry::operator!=(const Geometry& geometry) const
{
    return !(*this == geometry);
}

// ------------------------------------------------------------------------- //
// Member functions
// ------------------------------------------------------------------------- //

bool Geometry::IsInside(const Vector3D& position, const Vector3D& direction)
{
    bool is_inside = false;

    std::pair<double, double> dist = DistanceToBorder(position, direction);

    if (dist.first > 0 && dist.second < 0)
    {
        is_inside = true;
    }
    return is_inside;
}

// ------------------------------------------------------------------------- //
bool Geometry::IsInfront(const Vector3D& position, const Vector3D& direction)
{
    bool is_infront = false;

    std::pair<double, double> dist = DistanceToBorder(position, direction);

    if (dist.first > 0 && dist.second > 0)
    {
        is_infront = true;
    }
    return is_infront;
}

// ------------------------------------------------------------------------- //
bool Geometry::IsBehind(const Vector3D& position, const Vector3D& direction)
{
    bool is_behind = false;

    std::pair<double, double> dist = DistanceToBorder(position, direction);

    if (dist.first < 0 && dist.second < 0)
    {
        is_behind = true;
    }
    return is_behind;
}

// ------------------------------------------------------------------------- //
double Geometry::DistanceToClosestApproach(const Vector3D& position, const Vector3D& direction)
{
    return scalar_product(position_ - position, direction);
}
