
#include <algorithm>
#include <functional>
#include <cctype>
#include <string>
#include <memory>

#include "PROPOSAL/Logging.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"

namespace PROPOSAL {
static std::map<const std::string, std::shared_ptr<Medium>> Medium_Map
{
    {"water", std::shared_ptr<Medium>(new Water)},
    {"ice", std::shared_ptr<Medium>(new Ice)},
    {"salt", std::shared_ptr<Medium>(new Salt)},
    {"standardrock", std::shared_ptr<Medium>(new StandardRock)},
    {"frejusrock", std::shared_ptr<Medium>(new FrejusRock)},
    {"iron", std::shared_ptr<Medium>(new Iron)},
    {"hydrogen", std::shared_ptr<Medium>(new Hydrogen)},
    {"lead", std::shared_ptr<Medium>(new Lead)},
    {"copper", std::shared_ptr<Medium>(new Copper)},
    {"uranium", std::shared_ptr<Medium>(new Uranium)},
    {"air", std::shared_ptr<Medium>(new Air)},
    {"paraffin", std::shared_ptr<Medium>(new Paraffin)},
    {"antareswater", std::shared_ptr<Medium>(new AntaresWater)},
    {"cascadiabasinwater", std::shared_ptr<Medium>(new CascadiaBasinWater)},
};

std::shared_ptr<const Medium> CreateMedium(std::string name, double density_correction)
{
    std::transform(name.begin(), name.end(), name.begin(),
        [](unsigned char c){ return std::tolower(c); });

    std::unique_ptr<Density_distr> density_distr(new Density_homogeneous(density_correction));
    auto searched_medium = Medium_Map.find(name);
    if (searched_medium != Medium_Map.end()) {
        searched_medium->second->SetDensityDistribution(*density_distr);
    } else {
        log_fatal("Medium %s not registered!", name);
    }

    return searched_medium->second->create();
}
} // namespace PROPOSAL
