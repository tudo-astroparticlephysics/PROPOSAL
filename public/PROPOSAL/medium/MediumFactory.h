#include <memory>
#include <string>

#include "PROPOSAL/medium/Medium.h"

namespace PROPOSAL {
static std::map<const Medium_Type, std::shared_ptr<Medium>> Medium_Map{
    { Medium_Type::WATER, std::make_shared<Water>() },
    { Medium_Type::ICE, std::make_shared<Ice>() },
    { Medium_Type::SALT, std::make_shared<Salt>() },
    { Medium_Type::STANDARDROCK, std::make_shared<StandardRock>() },
    { Medium_Type::FREJUSROCK, std::make_shared<FrejusRock>() },
    { Medium_Type::IRON, std::make_shared<Iron>() },
    { Medium_Type::HYDROGEN, std::make_shared<Hydrogen>() },
    { Medium_Type::LEAD, std::make_shared<Lead>() },
    { Medium_Type::COPPER, std::make_shared<Copper>() },
    { Medium_Type::URANIUM, std::make_shared<Uranium>() },
    { Medium_Type::AIR, std::make_shared<Air>() },
    { Medium_Type::PARAFFIN, std::make_shared<Paraffin>() },
    { Medium_Type::ANTARESWATER, std::make_shared<AntaresWater>() },
    { Medium_Type::CASCADIABASINWATER, std::make_shared<CascadiaBasinWater>() },
};
} // namespace PROPOSAL

namespace PROPOSAL {
std::shared_ptr<Medium> CreateMedium(
    Medium_Type type, double density_correction = 1.0);
std::shared_ptr<Medium> CreateMedium(
    std::string name, double density_correction = 1.0);
} // namespace PROPOSAL
