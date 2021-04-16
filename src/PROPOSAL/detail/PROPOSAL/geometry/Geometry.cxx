/*
 * Geometry.cxx
 *
 *  Created on: 05.06.2013
 *      Author: koehne
 */

#include <sstream>
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

Geometry::Geometry(const std::string name, const Vector3D& position)
    : position_(position)
    , name_(name)
    , hierarchy_(0)
{
}

Geometry::Geometry(const nlohmann::json& config)
{
    if(!config.is_object()) throw std::invalid_argument("No json object found.");

    name_ = config.value("shape", "unknown");
    hierarchy_ = config.value("hierarchy", 0);

    if(!config.contains("origin"))
        throw std::invalid_argument("No geometry origin found.");
    position_ = Cartesian3D(config.at("origin"));
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

bool Geometry::IsInside(const Vector3D& position, const Vector3D& direction) const
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
bool Geometry::IsInfront(const Vector3D& position, const Vector3D& direction) const
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
bool Geometry::IsBehind(const Vector3D& position, const Vector3D& direction) const
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
bool Geometry::IsEntering(const Vector3D &position, const Vector3D &direction) const {
    auto dist_forward = DistanceToBorder(position, direction);
    auto dist_backward = DistanceToBorder(position, -Cartesian3D(direction));
    if (dist_forward.first >= 0 && dist_forward.second == -1) {
        if (dist_backward.first == -1 && dist_backward.second == -1) {
            return true;
        }
    }
    return false;
}

// ------------------------------------------------------------------------- //
bool Geometry::IsLeaving(const Vector3D &position, const Vector3D &direction) const {
    auto dist_forward = DistanceToBorder(position, direction);
    auto dist_backward = DistanceToBorder(position, -Cartesian3D(direction));
    if (dist_forward.first == -1 && dist_forward.second == -1) {
        if (dist_backward.first >= 0 && dist_backward.second == -1) {
            return true;
        }
    }
    return false;
}

Geometry::ParticleLocation::Enum Geometry::GetLocation(const Vector3D& position, const Vector3D& direction) const {
    if(IsInfront(position, direction))
        return Geometry::ParticleLocation::InfrontGeometry;
    if(IsInside(position, direction))
        return Geometry::ParticleLocation::InsideGeometry;
    else
        return Geometry::ParticleLocation::BehindGeometry;
}

// ------------------------------------------------------------------------- //
double Geometry::DistanceToClosestApproach(const Vector3D& position, const Vector3D& direction) const
{
    return (position_ - position) * direction;
}
