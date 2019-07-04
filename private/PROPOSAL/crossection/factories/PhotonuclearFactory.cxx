
#include <algorithm>

#include "PROPOSAL/crossection/PhotoIntegral.h"
#include "PROPOSAL/crossection/PhotoInterpolant.h"
#include "PROPOSAL/crossection/factories/PhotonuclearFactory.h"
#include "PROPOSAL/crossection/parametrization/PhotoQ2Integration.h"
#include "PROPOSAL/crossection/parametrization/PhotoRealPhotonAssumption.h"
#include "PROPOSAL/crossection/parametrization/Photonuclear.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

PhotonuclearFactory::PhotonuclearFactory()
    : photo_shadow_map_str_()
    , photo_shadow_map_enum_()
    , photo_real_map_str_()
    , photo_real_map_enum_()
    , photo_q2_map_str_()
    , photo_q2_map_enum_()
    , string_enum_()
    , string_shadow_enum_()
{
    // Register all photonuclear parametrizations in lower case!

    RegisterShadowEffect(
        "shadowduttarenosarcevicseckel", ShadowDuttaRenoSarcevicSeckel, &ShadowDuttaRenoSarcevicSeckel::create);
    RegisterShadowEffect("shadowbutkevichmikhailov", ShadowButkevichMikhailov, &ShadowButkevichMikhailov::create);

    RegisterRealPhoton("photozeus", Zeus, &PhotoZeus::create);
    RegisterRealPhoton("photobezrukovbugaev", BezrukovBugaev, &PhotoBezrukovBugaev::create);
    RegisterRealPhoton("photorhode", Rhode, &PhotoRhode::create);
    RegisterRealPhoton("photokokoulin", Kokoulin, &PhotoKokoulin::create);

    RegisterQ2("photoabramowiczlevinlevymaor91",
               AbramowiczLevinLevyMaor91,
               std::make_pair(&PhotoAbramowiczLevinLevyMaor91::create,
                              &PhotoQ2Interpolant<PhotoAbramowiczLevinLevyMaor91>::create));
    RegisterQ2("photoabramowiczlevinlevymaor97",
               AbramowiczLevinLevyMaor97,
               std::make_pair(&PhotoAbramowiczLevinLevyMaor97::create,
                              &PhotoQ2Interpolant<PhotoAbramowiczLevinLevyMaor97>::create));
    RegisterQ2("photobutkevichmikhailov",
               ButkevichMikhailov,
               std::make_pair(&PhotoButkevichMikhailov::create, &PhotoQ2Interpolant<PhotoButkevichMikhailov>::create));
    RegisterQ2("photorenosarcevicsu",
               RenoSarcevicSu,
               std::make_pair(&PhotoRenoSarcevicSu::create, &PhotoQ2Interpolant<PhotoRenoSarcevicSu>::create));
}

PhotonuclearFactory::~PhotonuclearFactory()
{
    photo_shadow_map_str_.clear();
    photo_shadow_map_enum_.clear();

    photo_real_map_str_.clear();
    photo_real_map_enum_.clear();

    photo_q2_map_str_.clear();
    photo_q2_map_enum_.clear();

    string_enum_.clear();
    string_shadow_enum_.clear();
}

// ------------------------------------------------------------------------- //
void PhotonuclearFactory::RegisterShadowEffect(const std::string& name,
                                               const Shadow& shadow,
                                               RegisterShadowEffectFunction create)
{
    photo_shadow_map_str_[name]    = create;
    photo_shadow_map_enum_[shadow] = create;
    string_shadow_enum_.insert(BimapStringShadowEnum::value_type(name, shadow));
}

// ------------------------------------------------------------------------- //
void PhotonuclearFactory::RegisterRealPhoton(const std::string& name, Enum enum_t, RegisterRealPhotonFunction create)
{
    photo_real_map_str_[name]    = create;
    photo_real_map_enum_[enum_t] = create;
    string_enum_.insert(BimapStringEnum::value_type(name, enum_t));
}

// ------------------------------------------------------------------------- //
void PhotonuclearFactory::RegisterQ2(const std::string& name,
                                     Enum enum_t,
                                     std::pair<RegisterQ2Function, RegisterQ2FunctionInterpolant> create)
{
    photo_q2_map_str_[name]    = create;
    photo_q2_map_enum_[enum_t] = create;
    string_enum_.insert(BimapStringEnum::value_type(name, enum_t));
}

// ------------------------------------------------------------------------- //
ShadowEffect* PhotonuclearFactory::CreateShadowEffect(const std::string& name)
{
    std::string name_lower = name;
    std::transform(name.begin(), name.end(), name_lower.begin(), ::tolower);

    PhotoShadowEffectMapString::const_iterator it = photo_shadow_map_str_.find(name_lower);

    if (it != photo_shadow_map_str_.end())
    {
        return it->second();
    } else
    {
        log_fatal("Photonuclear %s not registerd!", name.c_str());
        return NULL; // Just to prevent warinngs
    }
}

// ------------------------------------------------------------------------- //
ShadowEffect* PhotonuclearFactory::CreateShadowEffect(const Shadow& shadow)
{
    PhotoShadowEffectMapEnum::const_iterator it = photo_shadow_map_enum_.find(shadow);

    if (it != photo_shadow_map_enum_.end())
    {
        return it->second();
    } else
    {
        log_fatal("Photonuclear %s not registerd!", typeid(shadow).name());
        return NULL; // Just to prevent warinngs
    }
}

// ------------------------------------------------------------------------- //
CrossSection* PhotonuclearFactory::CreatePhotonuclear(const ParticleDef& particle_def,
                                                      const Medium& medium,
                                                      const EnergyCutSettings& cuts,
                                                      const Definition& def) const
{
    PhotoQ2MapEnum::const_iterator it_q2            = photo_q2_map_enum_.find(def.parametrization);
    PhotoRealPhotonMapEnum::const_iterator it_photo = photo_real_map_enum_.find(def.parametrization);

    if (it_q2 != photo_q2_map_enum_.end())
    {
        ShadowEffect* shadow = Get().CreateShadowEffect(def.shadow);

        PhotoIntegral* photo =
            new PhotoIntegral(*it_q2->second.first(particle_def, medium, cuts, def.multiplier, *shadow));
        delete shadow;
        return photo;
    } else if (it_photo != photo_real_map_enum_.end())
    {
        return new PhotoIntegral(*it_photo->second(particle_def, medium, cuts, def.multiplier, def.hard_component));
    } else
    {
        log_fatal("Photonuclear %s not registerd!", typeid(def.parametrization).name());
        return NULL; // Just to prevent warinngs
    }
}

// ------------------------------------------------------------------------- //
CrossSection* PhotonuclearFactory::CreatePhotonuclear(const ParticleDef& particle_def,
                                                      const Medium& medium,
                                                      const EnergyCutSettings& cuts,
                                                      const Definition& def,
                                                      InterpolationDef interpolation_def) const
{
    PhotoQ2MapEnum::const_iterator it_q2            = photo_q2_map_enum_.find(def.parametrization);
    PhotoRealPhotonMapEnum::const_iterator it_photo = photo_real_map_enum_.find(def.parametrization);

    if (it_q2 != photo_q2_map_enum_.end())
    {
        ShadowEffect* shadow = Get().CreateShadowEffect(def.shadow);

        PhotoInterpolant* photo = new PhotoInterpolant(
            *it_q2->second.second(particle_def, medium, cuts, def.multiplier, *shadow, interpolation_def),
            interpolation_def);
        delete shadow;
        return photo;
    } else if (it_photo != photo_real_map_enum_.end())
    {
        return new PhotoInterpolant(*it_photo->second(particle_def, medium, cuts, def.multiplier, def.hard_component),
                                    interpolation_def);
    } else
    {
        log_fatal("Photonuclear %s not registerd!", typeid(def.parametrization).name());
        return NULL; // Just to prevent warinngs
    }
}

// ------------------------------------------------------------------------- //
PhotonuclearFactory::Enum PhotonuclearFactory::GetEnumFromString(const std::string& name)
{
    std::string name_lower = name;
    std::transform(name.begin(), name.end(), name_lower.begin(), ::tolower);

    BimapStringEnum::left_const_iterator it = string_enum_.left.find(name_lower);
    if (it != string_enum_.left.end())
    {
        return it->second;
    } else
    {
        log_fatal("Photonuclear %s not registerd!", name.c_str());
        return PhotonuclearFactory::None; // Just to prevent warinngs
    }
}

// ------------------------------------------------------------------------- //
std::string PhotonuclearFactory::GetStringFromEnum(const PhotonuclearFactory::Enum& enum_t)
{
    BimapStringEnum::right_const_iterator it = string_enum_.right.find(enum_t);
    if (it != string_enum_.right.end())
    {
        return it->second;
    } else
    {
        log_fatal("Photonuclear %s not registerd!", typeid(enum_t).name());
        return ""; // Just to prevent warinngs
    }
}

// ------------------------------------------------------------------------- //
PhotonuclearFactory::Shadow PhotonuclearFactory::GetShadowEnumFromString(const std::string& name)
{
    std::string name_lower = name;
    std::transform(name.begin(), name.end(), name_lower.begin(), ::tolower);

    BimapStringShadowEnum::left_const_iterator it = string_shadow_enum_.left.find(name_lower);
    if (it != string_shadow_enum_.left.end())
    {
        return it->second;
    } else
    {
        log_fatal("Photonuclear %s not registerd!", name.c_str());
        return PhotonuclearFactory::ShadowNone; // Just to prevent warinngs
    }
}

// ------------------------------------------------------------------------- //
std::string PhotonuclearFactory::GetStringFromShadowEnum(const Shadow& shadow)
{

    BimapStringShadowEnum::right_const_iterator it = string_shadow_enum_.right.find(shadow);
    if (it != string_shadow_enum_.right.end())
    {
        return it->second;
    } else
    {
        log_fatal("Photonuclear %s not registerd!", typeid(shadow).name());
        return ""; // Just to prevent warinngs
    }
}
