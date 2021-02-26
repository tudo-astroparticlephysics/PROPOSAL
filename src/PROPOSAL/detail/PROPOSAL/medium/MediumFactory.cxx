
#include <algorithm>
#include <functional>
#include <cctype>
#include <string>
#include <memory>

#include "PROPOSAL/Logging.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"


namespace PROPOSAL {
std::shared_ptr<Medium> CreateMedium(Medium_Type type)
{
    auto searched_medium = Medium_Map.find(type);
    if (searched_medium == Medium_Map.end()) {
        throw std::invalid_argument("Medium not found.");
    }
    return searched_medium->second;
}
} // namespace PROPOSAL

namespace PROPOSAL {
std::shared_ptr<Medium> CreateMedium(std::string name)
{
    std::transform(name.begin(), name.end(), name.begin(),
        [](unsigned char c){ return std::tolower(c); });

    for (size_t id = 0; id < Medium_Name.size(); ++id) {
        if (name ==Medium_Name[id]){
            return CreateMedium(static_cast<Medium_Type>(id));
        }
    }
    throw std::invalid_argument("Medium not found.");
}
} // namespace PROPOSAL