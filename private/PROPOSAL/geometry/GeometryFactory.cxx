
#include <algorithm>
#include <memory>
#include <stdexcept>

#include "PROPOSAL/Logging.h"
#include "PROPOSAL/geometry/Geometry.h"
#include "PROPOSAL/geometry/GeometryFactory.h"

namespace PROPOSAL {
std::shared_ptr<Geometry> CreateGeometry(Geometry_Type type)
{
    auto searched_geometry = Geometry_Map.find(type);
    if (searched_geometry != Geometry_Map.end()) {
        return searched_geometry->second;
    }
    throw std::invalid_argument("Geometry not found.");
}
} // namespace PROPOSAL

namespace PROPOSAL {
std::shared_ptr<Geometry> CreateGeometry(std::string name)
{
    std::transform(name.begin(), name.end(), name.begin(),
        [](unsigned char c) { return std::tolower(c); });

    for (size_t id = 0; id < Geometry_Name.size(); ++id) {
        if (name == Geometry_Name[id]) {
            return CreateGeometry(static_cast<Geometry_Type>(id));
        }
    }
    throw std::invalid_argument("Geometry not found.");
}
} // namespace PROPOSAL
