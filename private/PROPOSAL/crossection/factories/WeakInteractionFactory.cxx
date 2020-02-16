
#include <algorithm>

#include "PROPOSAL/crossection/WeakIntegral.h"
#include "PROPOSAL/crossection/WeakInterpolant.h"
#include "PROPOSAL/crossection/factories/WeakInteractionFactory.h"
#include "PROPOSAL/crossection/parametrization/WeakInteraction.h"

#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

WeakInteractionFactory::WeakInteractionFactory()
        : weak_map_str_()
        , weak_map_enum_()
        , string_enum_()
{
    Register("weakcoopersarkarmertsch", CooperSarkarMertsch, &WeakCooperSarkarMertsch::create);
    Register("none", None, nullptr);
}

WeakInteractionFactory::~WeakInteractionFactory()
{
    string_enum_.clear();
    weak_map_str_.clear();
    weak_map_enum_.clear();
}

// ------------------------------------------------------------------------- //
// Most general creator
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
CrossSection* WeakInteractionFactory::CreateWeakInteraction(const ParticleDef& particle_def,
                                                            std::shared_ptr<const Medium> medium,
                                                            const Definition& def) const
{
    if(def.parametrization == WeakInteractionFactory::Enum::None){
        log_fatal("Can't return Weakinteraction Crosssection if parametrization is None");
        return NULL;
    }

    WeakMapEnum::const_iterator it = weak_map_enum_.find(def.parametrization);

    if (it != weak_map_enum_.end())
    {
        return new WeakIntegral(*it->second(particle_def, medium, def.multiplier));
    } else
    {
        log_fatal("WeakInteraction %s not registered!", typeid(def.parametrization).name());
        return NULL; // Just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
CrossSection* WeakInteractionFactory::CreateWeakInteraction(const ParticleDef& particle_def,
                                                            std::shared_ptr<const Medium> medium,
                                                            const Definition& def,
                                                            InterpolationDef interpolation_def) const
{
    if(def.parametrization == WeakInteractionFactory::Enum::None){
        log_fatal("Can't return Weakinteraction Crosssection if parametrization is None");
        return NULL;
    }

    WeakMapEnum::const_iterator it = weak_map_enum_.find(def.parametrization);

    if (it != weak_map_enum_.end())
    {
        return new WeakInterpolant(*it->second(particle_def, medium, def.multiplier), interpolation_def);
    } else
    {
        log_fatal("WeakInteraction %s not registered!", typeid(def.parametrization).name());
        return NULL; // Just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
void WeakInteractionFactory::Register(const std::string& name,
                                      Enum enum_t,
                                      RegisterFunction create)
{
    weak_map_str_[name]    = create;
    weak_map_enum_[enum_t] = create;
    string_enum_.insert(name, enum_t);
}

// ------------------------------------------------------------------------- //
WeakInteractionFactory::Enum WeakInteractionFactory::GetEnumFromString(const std::string& name)
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
        log_fatal("WeakInteraction %s not registered!", name.c_str());
        return Fail; // Just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
std::string WeakInteractionFactory::GetStringFromEnum(const WeakInteractionFactory::Enum& enum_t)
{
    auto& right = string_enum_.GetRight();
    auto it = right.find(enum_t);
    if (it != right.end())
    {
        return it->second;
    } else
    {
        log_fatal("WeakInteraction %s not registered!", typeid(enum_t).name());
        return ""; // Just to prevent warnings
    }
}
