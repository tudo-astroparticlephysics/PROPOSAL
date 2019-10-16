
#include <algorithm>

#include "PROPOSAL/crossection/EpairIntegral.h"
#include "PROPOSAL/crossection/EpairInterpolant.h"
#include "PROPOSAL/crossection/factories/EpairProductionFactory.h"
#include "PROPOSAL/crossection/parametrization/EpairProduction.h"

#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

EpairProductionFactory::EpairProductionFactory()
    : epair_map_str_()
    , epair_map_enum_()
    , string_enum_()
{
    Register("epairkelnerkokoulinpetrukhin", KelnerKokoulinPetrukhin, std::make_pair(&EpairKelnerKokoulinPetrukhin::create, &EpairProductionRhoInterpolant<EpairKelnerKokoulinPetrukhin>::create));
    Register("epairsandrocksoedingreksorhode", SandrockSoedingreksoRhode, std::make_pair(&EpairSandrockSoedingreksoRhode::create, &EpairProductionRhoInterpolant<EpairSandrockSoedingreksoRhode>::create));
}

EpairProductionFactory::~EpairProductionFactory()
{
    string_enum_.clear();
    epair_map_str_.clear();
    epair_map_enum_.clear();
}

// ------------------------------------------------------------------------- //
// Most general creator
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
CrossSection* EpairProductionFactory::CreateEpairProduction(const ParticleDef& particle_def,
                                                            const Medium& medium,
                                                            const EnergyCutSettings& cuts,
                                                            const Definition& def) const
{
    EpairMapEnum::const_iterator it = epair_map_enum_.find(def.parametrization);

    if (it != epair_map_enum_.end())
    {
        return new EpairIntegral(*it->second.first(particle_def, medium, cuts, def.multiplier, def.lpm_effect));
    } else
    {
        log_fatal("EpairProduction %s not registerd!", typeid(def.parametrization).name());
        return NULL; // Just to prevent warinngs
    }
}

// ------------------------------------------------------------------------- //
CrossSection* EpairProductionFactory::CreateEpairProduction(const ParticleDef& particle_def,
                                                            const Medium& medium,
                                                            const EnergyCutSettings& cuts,
                                                            const Definition& def,
                                                            InterpolationDef interpolation_def) const
{
    EpairMapEnum::const_iterator it = epair_map_enum_.find(def.parametrization);

    if (it != epair_map_enum_.end())
    {
        return new EpairInterpolant(*it->second.second(particle_def, medium, cuts, def.multiplier, def.lpm_effect, interpolation_def), interpolation_def);
    } else
    {
        log_fatal("EpairProduction %s not registerd!", typeid(def.parametrization).name());
        return NULL; // Just to prevent warinngs
    }
}

// ------------------------------------------------------------------------- //
void EpairProductionFactory::Register(const std::string& name,
                                     Enum enum_t,
                                     std::pair<RegisterFunction, RegisterFunctionInterpolant> create)
{
    epair_map_str_[name]    = create;
    epair_map_enum_[enum_t] = create;
    string_enum_.insert(name, enum_t);
}

// ------------------------------------------------------------------------- //
EpairProductionFactory::Enum EpairProductionFactory::GetEnumFromString(const std::string& name)
{
    std::string name_lower = name;
    std::transform(name.begin(), name.end(), name_lower.begin(), ::tolower);

    auto& left = string_enum_.GetLeft();
    auto it = left.find(name_lower);
    if (it != left.end())
    {
        return it->second;
    } else
    {
        log_fatal("EpairProduction %s not registerd!", name.c_str());
        return None; // Just to prevent warinngs
    }
}

// ------------------------------------------------------------------------- //
std::string EpairProductionFactory::GetStringFromEnum(const EpairProductionFactory::Enum& enum_t)
{
    auto& right = string_enum_.GetRight();
    auto it = right.find(enum_t);
    if (it != right.end())
    {
        return it->second;
    } else
    {
        log_fatal("EpairProduction %s not registerd!", typeid(enum_t).name());
        return ""; // Just to prevent warinngs
    }
}
