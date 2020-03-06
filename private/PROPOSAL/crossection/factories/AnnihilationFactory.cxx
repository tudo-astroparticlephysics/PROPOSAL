
#include <algorithm>

#include "PROPOSAL/crossection/AnnihilationIntegral.h"
#include "PROPOSAL/crossection/AnnihilationInterpolant.h"
#include "PROPOSAL/crossection/factories/AnnihilationFactory.h"
#include "PROPOSAL/crossection/parametrization/Annihilation.h"

#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

AnnihilationFactory::AnnihilationFactory()
        : annihilation_map_str_()
        , annihilation_map_enum_()
        , string_enum_()
{
    Register("annihilationheitler", Heitler, &AnnihilationHeitler::create);
    Register("none", None, nullptr);
}

AnnihilationFactory::~AnnihilationFactory()
{
    string_enum_.clear();
    annihilation_map_str_.clear();
    annihilation_map_enum_.clear();
}

// ------------------------------------------------------------------------- //
// Most general creator
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
CrossSection* AnnihilationFactory::CreateAnnihilation(const ParticleDef& particle_def,
                                                            std::shared_ptr<const Medium> medium,
                                                            const Definition& def) const
{
    if(def.parametrization == AnnihilationFactory::Enum::None){
        log_fatal("Can't return Annihilation Crosssection if parametrization is None");
        return NULL;
    }

    AnnihilationMapEnum::const_iterator it = annihilation_map_enum_.find(def.parametrization);

    if (it != annihilation_map_enum_.end())
    {
        return new AnnihilationIntegral(*it->second(particle_def, medium, def.multiplier));
    } else
    {
        log_fatal("Annihilation %s not registered!", typeid(def.parametrization).name());
        return NULL; // Just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
CrossSection* AnnihilationFactory::CreateAnnihilation(const ParticleDef& particle_def,
                                                            std::shared_ptr<const Medium> medium,
                                                            const Definition& def,
                                                            InterpolationDef interpolation_def) const
{
    if(def.parametrization == AnnihilationFactory::Enum::None){
        log_fatal("Can't return Annihilation Crosssection if parametrization is None");
        return NULL;
    }

    AnnihilationMapEnum::const_iterator it = annihilation_map_enum_.find(def.parametrization);

    if (it != annihilation_map_enum_.end())
    {
        return new AnnihilationInterpolant(*it->second(particle_def, medium, def.multiplier), interpolation_def);
    } else
    {
        log_fatal("Annihilation %s not registered!", typeid(def.parametrization).name());
        return NULL; // Just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
void AnnihilationFactory::Register(const std::string& name,
                                      Enum enum_t,
                                      RegisterFunction create)
{
    annihilation_map_str_[name]    = create;
    annihilation_map_enum_[enum_t] = create;
    string_enum_.insert(name, enum_t);
}

// ------------------------------------------------------------------------- //
AnnihilationFactory::Enum AnnihilationFactory::GetEnumFromString(const std::string& name)
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
        log_fatal("Annihilation %s not registered!", name.c_str());
        return AnnihilationFactory::Fail; // Just to prevent warnings

    }
}

// ------------------------------------------------------------------------- //
std::string AnnihilationFactory::GetStringFromEnum(const AnnihilationFactory::Enum& enum_t)
{
    auto& right = string_enum_.GetRight();
    auto it = right.find(enum_t);
    if (it != right.end())
    {
        return it->second;
    } else
    {
        log_fatal("Annihilation %s not registered!", typeid(enum_t).name());
        return ""; // Just to prevent warnings
    }
}
