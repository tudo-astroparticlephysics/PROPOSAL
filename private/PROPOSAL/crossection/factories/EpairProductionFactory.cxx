
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
    Register("epairkelnerkokoulinpetrukhin", KelnerKokoulinPetrukhin, &EpairKelnerKokoulinPetrukhin::create);
    Register("epairsandrocksoedingreksorhode", SandrockSoedingreksoRhode, &EpairSandrockSoedingreksoRhode::create);
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
                                                            std::shared_ptr<const Medium> medium,
                                                            std::shared_ptr<const EnergyCutSettings> cuts,
                                                            const Definition& def,
                                                            std::shared_ptr<const InterpolationDef> interpolation_def) const
{
    EpairMapEnum::const_iterator it = epair_map_enum_.find(def.parametrization);

    if (it != epair_map_enum_.end())
    {
        if(interpolation_def == nullptr){
            return new EpairIntegral(*it->second(particle_def, medium, def.multiplier, def.lpm_effect), cuts);
        }
        else{
            return new EpairInterpolant(*it->second(particle_def, medium, def.multiplier, def.lpm_effect), cuts, *interpolation_def);
        }
    } else
    {
        log_fatal("EpairProduction %s not registered!", typeid(def.parametrization).name());
        return NULL; // Just to prevent warnings
    }
}

CrossSection* EpairProductionFactory::CreateEpairProduction(const EpairProduction& parametrization,
                                                            std::shared_ptr<const EnergyCutSettings> cuts,
                                                            std::shared_ptr<const InterpolationDef> interpolation_def) const
{
    if(interpolation_def==nullptr){
        return new EpairIntegral(parametrization, cuts);
    }
    else{
        return new EpairInterpolant(parametrization, cuts, *interpolation_def);
    }
}

// ------------------------------------------------------------------------- //
void EpairProductionFactory::Register(const std::string& name,
                                     Enum enum_t,
                                     RegisterFunction create)
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
        log_fatal("EpairProduction %s not registered!", name.c_str());
        return Fail; // Just to prevent warnings
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
        log_fatal("EpairProduction %s not registered!", typeid(enum_t).name());
        return ""; // Just to prevent warnings
    }
}
