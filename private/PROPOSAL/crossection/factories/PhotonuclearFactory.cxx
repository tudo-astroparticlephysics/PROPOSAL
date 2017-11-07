
#include <boost/algorithm/string.hpp> // case insensitive string compare for configuration file

#include "PROPOSAL/crossection/PhotoIntegral.h"
#include "PROPOSAL/crossection/PhotoInterpolant.h"
#include "PROPOSAL/crossection/parametrization/Photonuclear.h"
#include "PROPOSAL/crossection/parametrization/PhotoRealPhotonAssumption.h"
#include "PROPOSAL/crossection/parametrization/PhotoQ2Integration.h"
#include "PROPOSAL/crossection/factories/PhotonuclearFactory.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/Output.h"

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

    RegisterShadowEffect("drss", ShadowDuttaRenoSarcevicSeckel, &ShadowDuttaRenoSarcevicSeckel::create);
    RegisterShadowEffect("butkevich-mikhailov", ShadowButkevichMikhailov, &ShadowButkevichMikhailov::create);

    RegisterRealPhoton("zeus", Zeus, &PhotoZeus::create);
    RegisterRealPhoton("bezrukov-bugaev", BezrukovBugaev, &PhotoBezrukovBugaev::create);
    RegisterRealPhoton("rhode", Rhode, &PhotoRhode::create);
    RegisterRealPhoton("kokoulin", Kokoulin, &PhotoKokoulin::create);

    RegisterQ2("allm91",
             AbramowiczLevinLevyMaor91,
             std::make_pair(&PhotoAbramowiczLevinLevyMaor91::create,
                            &PhotoQ2Interpolant<PhotoAbramowiczLevinLevyMaor91>::create));
    RegisterQ2("allm97",
             AbramowiczLevinLevyMaor97,
             std::make_pair(&PhotoAbramowiczLevinLevyMaor97::create,
                            &PhotoQ2Interpolant<PhotoAbramowiczLevinLevyMaor97>::create));
    RegisterQ2("butkevich-mikhailov",
             ButkevichMikhailov,
             std::make_pair(&PhotoButkevichMikhailov::create,
                            &PhotoQ2Interpolant<PhotoButkevichMikhailov>::create));
    RegisterQ2("reno-sarcevic-su",
             RenoSarcevicSu,
             std::make_pair(&PhotoRenoSarcevicSu::create,
                            &PhotoQ2Interpolant<PhotoRenoSarcevicSu>::create));
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
void PhotonuclearFactory::RegisterShadowEffect(const std::string& name, const Shadow& shadow, RegisterShadowEffectFunction create)
{
    photo_shadow_map_str_[name] = create;
    photo_shadow_map_enum_[shadow] = create;
    string_shadow_enum_.insert(BimapStringShadowEnum::value_type(name, shadow));
}

// ------------------------------------------------------------------------- //
void PhotonuclearFactory::RegisterRealPhoton(const std::string& name, Enum enum_t, RegisterRealPhotonFunction create)
{
    photo_real_map_str_[name] = create;
    photo_real_map_enum_[enum_t] = create;
    string_enum_.insert(BimapStringEnum::value_type(name, enum_t));
}

// ------------------------------------------------------------------------- //
void PhotonuclearFactory::RegisterQ2(const std::string& name, Enum enum_t, std::pair<RegisterQ2Function, RegisterQ2FunctionInterpolant> create)
{
    photo_q2_map_str_[name] = create;
    photo_q2_map_enum_[enum_t] = create;
    string_enum_.insert(BimapStringEnum::value_type(name, enum_t));
}

// ------------------------------------------------------------------------- //
ShadowEffect* PhotonuclearFactory::CreateShadowEffect(const std::string& name)
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    PhotoShadowEffectMapString::const_iterator it = photo_shadow_map_str_.find(name_lower);

    if (it != photo_shadow_map_str_.end())
    {
        return it->second();
    } else
    {
        log_fatal("Photonuclear %s not registerd!", name.c_str());
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
    }
}

// ------------------------------------------------------------------------- //
CrossSection* PhotonuclearFactory::CreatePhotonuclear(const ParticleDef& particle_def,
                                                      const Medium& medium,
                                                      const EnergyCutSettings& cuts,
                                                      const Definition& def) const
{
    PhotoQ2MapEnum::const_iterator it_q2 = photo_q2_map_enum_.find(def.parametrization);
    PhotoRealPhotonMapEnum::const_iterator it_photo = photo_real_map_enum_.find(def.parametrization);

    if (it_q2 != photo_q2_map_enum_.end())
    {
        ShadowEffect* shadow = Get().CreateShadowEffect(def.shadow);

        PhotoIntegral* photo =  new PhotoIntegral(*it_q2->second.first(particle_def, medium, cuts, *shadow, def.multiplier));
        delete shadow;
        return photo;
    }
    else if (it_photo != photo_real_map_enum_.end())
    {
        RealPhoton* real_photon;

        if (def.hardbb)
        {
            real_photon = new HardBB(particle_def);
        }
        else
        {
            real_photon = new SoftBB();
        }

        PhotoIntegral* photo = new PhotoIntegral(*it_photo->second(particle_def, medium, cuts, *real_photon, def.multiplier));
        delete real_photon;
        return photo;
    }
    else
    {
        log_fatal("Photonuclear %s not registerd!", typeid(def.parametrization).name());
    }
}

// ------------------------------------------------------------------------- //
CrossSection* PhotonuclearFactory::CreatePhotonuclear(const ParticleDef& particle_def,
                                                      const Medium& medium,
                                                      const EnergyCutSettings& cuts,
                                                      const Definition& def,
                                                      InterpolationDef interpolation_def) const
{
    PhotoQ2MapEnum::const_iterator it_q2 = photo_q2_map_enum_.find(def.parametrization);
    PhotoRealPhotonMapEnum::const_iterator it_photo = photo_real_map_enum_.find(def.parametrization);

    if (it_q2 != photo_q2_map_enum_.end())
    {
        ShadowEffect* shadow = Get().CreateShadowEffect(def.shadow);

        PhotoInterpolant* photo =  new PhotoInterpolant(*it_q2->second.second(particle_def, medium, cuts, *shadow, def.multiplier, interpolation_def), interpolation_def);
        delete shadow;
        return photo;
    }
    else if (it_photo != photo_real_map_enum_.end())
    {
        RealPhoton* real_photon;

        if (def.hardbb)
        {
            real_photon = new HardBB(particle_def);
        }
        else
        {
            real_photon = new SoftBB();
        }

        PhotoInterpolant* photo = new PhotoInterpolant(*it_photo->second(particle_def, medium, cuts, *real_photon, def.multiplier), interpolation_def);
        delete real_photon;
        return photo;
    }
    else
    {
        log_fatal("Photonuclear %s not registerd!", typeid(def.parametrization).name());
    }
}

// ------------------------------------------------------------------------- //
PhotonuclearFactory::Enum PhotonuclearFactory::GetEnumFromString(const std::string& name)
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    BimapStringEnum::left_const_iterator it = string_enum_.left.find(name_lower);
    if (it != string_enum_.left.end())
    {
        return it->second;
    } else
    {
        log_fatal("Photonuclear %s not registerd!", name.c_str());
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
    }
}

// ------------------------------------------------------------------------- //
PhotonuclearFactory::Shadow PhotonuclearFactory::GetShadowEnumFromString(const std::string& name)
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    BimapStringShadowEnum::left_const_iterator it = string_shadow_enum_.left.find(name_lower);
    if (it != string_shadow_enum_.left.end())
    {
        return it->second;
    } else
    {
        log_fatal("Photonuclear %s not registerd!", name.c_str());
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
    }
}
