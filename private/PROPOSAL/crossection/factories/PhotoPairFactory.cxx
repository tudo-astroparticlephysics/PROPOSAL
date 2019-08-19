
#include <algorithm>

#include "PROPOSAL/crossection/PhotoPairIntegral.h"
#include "PROPOSAL/crossection/PhotoPairInterpolant.h"
#include "PROPOSAL/crossection/factories/PhotoPairFactory.h"
#include "PROPOSAL/crossection/parametrization/PhotoPairProduction.h"

#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

PhotoPairFactory::PhotoPairFactory()
        : photopair_map_str_()
        , photopair_map_enum_()
        , string_enum_()
{
    // Register all PhotoPair parametrizations in lower case!

    Register("photopairtsai", Tsai, &PhotoPairTsai::create);
    Register("none", None, nullptr); //empty parametrization
}

// ------------------------------------------------------------------------- //
PhotoPairFactory::~PhotoPairFactory()
{
    photopair_map_str_.clear();
    photopair_map_enum_.clear();
    string_enum_.clear();
}

// ------------------------------------------------------------------------- //
void PhotoPairFactory::Register(const std::string& name, Enum enum_t, RegisterFunction create)
{
    photopair_map_str_[name]    = create;
    photopair_map_enum_[enum_t] = create;
    string_enum_.insert(BimapStringEnum::value_type(name, enum_t));
}

// --------------------------------------------------------------------- //
// Most general creation
// --------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
CrossSection* PhotoPairFactory::CreatePhotoPair(const ParticleDef& particle_def,
                                            const Medium& medium,
                                            const Definition& def) const
{
    if(def.parametrization == PhotoPairFactory::Enum::None){
        log_fatal("Can't return PhotoPair Crosssection if parametrization is None");
        return NULL;
    }

    PhotoPairMapEnum::const_iterator it = photopair_map_enum_.find(def.parametrization);

    if (it != photopair_map_enum_.end())
    {
        return new PhotoPairIntegral(*it->second(particle_def, medium, def.multiplier));
    } else
    {
        log_fatal("PhotoPair %s not registered!", typeid(def.parametrization).name());
        return NULL; // Just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
CrossSection* PhotoPairFactory::CreatePhotoPair(const ParticleDef& particle_def,
                                            const Medium& medium,
                                            const Definition& def,
                                            InterpolationDef interpolation_def) const
{
    if(def.parametrization == PhotoPairFactory::Enum::None){
        log_fatal("Can't return PhotoPair Crosssection if parametrization is None");
        return NULL;
    }

    PhotoPairMapEnum::const_iterator it = photopair_map_enum_.find(def.parametrization);

    if (it != photopair_map_enum_.end())
    {
        return new PhotoPairInterpolant(*it->second(particle_def, medium, def.multiplier), interpolation_def);
    } else
    {
        log_fatal("PhotoPair %s not registered!", typeid(def.parametrization).name());
        return NULL; // Just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
PhotoPairFactory::Enum PhotoPairFactory::GetEnumFromString(const std::string& name)
{
    std::string name_lower = name;
    std::transform(name.begin(), name.end(), name_lower.begin(), ::tolower);

    BimapStringEnum::left_const_iterator it = string_enum_.left.find(name_lower);
    if (it != string_enum_.left.end())
    {
        return it->second;
    } else
    {
        log_fatal("PhotoPair %s not registered!", name.c_str());
        return PhotoPairFactory::Fail; // Just to prevent warinngs

    }
}

// ------------------------------------------------------------------------- //
std::string PhotoPairFactory::GetStringFromEnum(const PhotoPairFactory::Enum& enum_t)
{
    BimapStringEnum::right_const_iterator it = string_enum_.right.find(enum_t);
    if (it != string_enum_.right.end())
    {
        return it->second;
    } else
    {
        log_fatal("PhotoPair %s not registered!", typeid(enum_t).name());
        return ""; // Just to prevent warnings
    }
}
