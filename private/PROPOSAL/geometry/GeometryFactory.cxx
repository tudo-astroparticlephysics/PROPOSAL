
#include <algorithm>
#include <memory>
#include <stdexcept>

#include "PROPOSAL/Logging.h"
#include "PROPOSAL/geometry/Geometry.h"
#include "PROPOSAL/geometry/GeometryFactory.h"

namespace PROPOSAL {
std::shared_ptr<Geometry> GetGeometry(std::string name)
{
    std::transform(name.begin(), name.end(), name.begin(),
        [](unsigned char c) { return std::tolower(c); });

    std::unique_ptr<Geometry> geometry;
    auto searched_geometry = Geometry_Map.find(name);
    if (searched_geometry != Geometry_Map.end()) {
        return searched_geometry->second;
    }

    throw std::invalid_argument("Geometry not found.");
}
} // namespace PROPOSAL

namespace PROPOSAL {
std::shared_ptr<const Geometry> CreateGeometry(std::string name)
{
    return GetGeometry(name)->create();
}
} // namespace PROPOSAL
