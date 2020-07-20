
#include <algorithm>
#include <functional>
#include <cctype>
#include <string>
#include <memory>

#include "PROPOSAL/Logging.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"


namespace PROPOSAL {
std::shared_ptr<Medium> CreateMedium(Medium_Type type, double density_correction)
{
    //std::unique_ptr<Density_distr> density_distr(new Density_homogeneous(density_correction));
    auto searched_medium = Medium_Map.find(type);
    if (searched_medium == Medium_Map.end()) {
        throw std::invalid_argument("Medium not found.");
    }
    throw std::invalid_argument("Medium not found.");
}
} // namespace PROPOSAL

namespace PROPOSAL {
std::shared_ptr<Medium> CreateMedium(std::string name, double density_correction)
{
    std::transform(name.begin(), name.end(), name.begin(),
        [](unsigned char c){ return std::tolower(c); });

    for (size_t id = 0; id < Medium_Name.size(); ++id) {
        if (name ==Medium_Name[id]){
            return CreateMedium(static_cast<Medium_Type>(id), density_correction);
        }
    }
    throw std::invalid_argument("Medium not found.");
}
} // namespace PROPOSAL

namespace PROPOSAL{
std::shared_ptr<Medium> CreateMedium(const nlohmann::json& config){
    if (config.contains("name"))
    {
        std::string name;
        config.at("name").get_to(name);
        auto density_correction = config.value("density_correction", 1.0);
        return CreateMedium(name, density_correction);
    }
    else {
        throw std::invalid_argument("Medium must be specified via parameter 'name'");
    }
}
} // namespace PROPOSAL
