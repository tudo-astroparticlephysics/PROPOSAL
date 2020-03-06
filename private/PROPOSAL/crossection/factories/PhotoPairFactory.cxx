
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
        , photo_angle_map_str_()
        , photo_angle_map_enum_()
        , string_photo_angle_enum_()
{
    // Register all PhotoPair parametrizations in lower case!

    Register("photopairtsai", Tsai, &PhotoPairTsai::create);
    Register("none", None, nullptr); //empty parametrization

    // Register PhotoAngle parametrizations

    RegisterPhotoAngle("photoangletsaiintegral", PhotoAngleTsaiIntegral, &PhotoAngleTsaiIntegral::create);
    RegisterPhotoAngle("photoanglenodeflection", PhotoAngleNoDeflection, &PhotoAngleNoDeflection::create);
    RegisterPhotoAngle("photoangleegs", PhotoAngleEGS, &PhotoAngleEGS::create);
}

// ------------------------------------------------------------------------- //
PhotoPairFactory::~PhotoPairFactory()
{
    photopair_map_str_.clear();
    photopair_map_enum_.clear();
    string_enum_.clear();

    //PhotoAngle
    photo_angle_map_str_.clear();
    photo_angle_map_enum_.clear();
    string_photo_angle_enum_.clear();
}

// ------------------------------------------------------------------------- //
void PhotoPairFactory::Register(const std::string& name, Enum enum_t, RegisterFunction create)
{
    photopair_map_str_[name]    = create;
    photopair_map_enum_[enum_t] = create;
    string_enum_.insert(name, enum_t);
}

// --------------------------------------------------------------------- //
// Most general creation
// --------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
CrossSection* PhotoPairFactory::CreatePhotoPair(const ParticleDef& particle_def,
                                            std::shared_ptr<const Medium> medium,
                                            const Definition& def) const
{
    if(def.parametrization == PhotoPairFactory::Enum::None){
        log_fatal("Can't return PhotoPair Crosssection if parametrization is None");
        return NULL;
    }

    PhotoPairMapEnum::const_iterator it = photopair_map_enum_.find(def.parametrization);

    if (it != photopair_map_enum_.end())
    {
        PhotoAngleDistribution* photoangle = Get().CreatePhotoAngleDistribution(def.photoangle, particle_def, medium);
        PhotoPairIntegral* photopair = new PhotoPairIntegral(*it->second(particle_def, medium, def.multiplier), *photoangle);
        delete photoangle;
        return photopair;
    } else
    {
        log_fatal("PhotoPair %s not registered!", typeid(def.parametrization).name());
        return NULL; // Just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
CrossSection* PhotoPairFactory::CreatePhotoPair(const ParticleDef& particle_def,
                                            std::shared_ptr<const Medium> medium,
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
        PhotoAngleDistribution* photoangle = Get().CreatePhotoAngleDistribution(def.photoangle, particle_def, medium);
        PhotoPairInterpolant* photopair = new PhotoPairInterpolant(*it->second(particle_def, medium, def.multiplier), *photoangle, interpolation_def);
        delete photoangle;
        return photopair;
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

    auto& left = string_enum_.GetLeft();
    auto it = left.find(name_lower);
    if (it != left.end())
    {
        return it->second;
    } else
    {
        log_fatal("PhotoPair %s not registered!", name.c_str());
        return PhotoPairFactory::Fail; // Just to prevent warnings

    }
}

// ------------------------------------------------------------------------- //
std::string PhotoPairFactory::GetStringFromEnum(const PhotoPairFactory::Enum& enum_t)
{
    auto& right = string_enum_.GetRight();
    auto it = right.find(enum_t);
    if (it != right.end())
    {
        return it->second;
    } else
    {
        log_fatal("PhotoPair %s not registered!", typeid(enum_t).name());
        return ""; // Just to prevent warnings
    }
}

// --------------------------------------------------------------------- //
// PhotoAngleDistribution
// --------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
void PhotoPairFactory::RegisterPhotoAngle(const std::string& name, const PhotoAngle& photoangle_t, RegisterPhotoAngleFunction create)
{
    photo_angle_map_str_[name]    = create;
    photo_angle_map_enum_[photoangle_t] = create;
    string_photo_angle_enum_.insert(name, photoangle_t);}

// ------------------------------------------------------------------------- //
PhotoAngleDistribution* PhotoPairFactory::CreatePhotoAngleDistribution(
        const std::string& name,
        const ParticleDef& particle_def,
        std::shared_ptr<const Medium> medium) {
    std::string name_lower = name;
    std::transform(name.begin(), name.end(), name_lower.begin(), ::tolower);

    PhotoAngleMapString::const_iterator it = photo_angle_map_str_.find(name_lower);

    if (it != photo_angle_map_str_.end())
    {
        return it->second(particle_def, medium);
    } else
    {
        log_fatal("PhotoAngleDistribution %s not registered!", name.c_str());
        return NULL; // Just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
PhotoAngleDistribution* PhotoPairFactory::CreatePhotoAngleDistribution(
        const PhotoPairFactory::PhotoAngle& photoangle,
        const ParticleDef& particle_def,
        std::shared_ptr<const Medium> medium) {

    PhotoAngleMapEnum::const_iterator it = photo_angle_map_enum_.find(photoangle);

    if (it != photo_angle_map_enum_.end())
    {
        return it->second(particle_def, medium);
    } else
    {
        log_fatal("PhotoAngleDistribution %s not registered!", typeid(photoangle).name());
        return NULL; // Just to prevent warnings
    }
}

// ------------------------------------------------------------------------- //
PhotoPairFactory::PhotoAngle PhotoPairFactory::GetPhotoAngleEnumFromString(const std::string& name)
{
    std::string name_lower = name;
    std::transform(name.begin(), name.end(), name_lower.begin(), ::tolower);

    auto& left = string_photo_angle_enum_.GetLeft();
    auto it = left.find(name_lower);
    if (it != left.end())
    {
        return it->second;
    } else
    {
        log_fatal("PhotoAngle %s not registered!", name.c_str());
        return PhotoPairFactory::PhotoAngleNoDeflection; // Just to prevent warnings

    }
}

// ------------------------------------------------------------------------- //
std::string PhotoPairFactory::GetStringFromPhotoAngleEnum(const PhotoAngle& photoangle)
{
    auto& right = string_photo_angle_enum_.GetRight();
    auto it = right.find(photoangle);
    if (it != right.end())
    {
        return it->second;
    } else
    {
        log_fatal("PhotoAngle %s not registered!", typeid(photoangle).name());
        return ""; // Just to prevent warnings
    }
}
