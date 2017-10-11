
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
{
    // Register all photonuclear parametrizations in lower case!

    RegisterShadowEffect("dutta", ShadowDutta, &ShadowDutta::create);
    RegisterShadowEffect("butkevichmikhailov", ShadowButkevichMikhailov, &ShadowButkevichMikhailov::create);

    RegisterRealPhoton("zeus", Zeus, &PhotoZeus::create);
    RegisterRealPhoton("bezrukovbugaev", BezrukovBugaev, &PhotoBezrukovBugaev::create);
    RegisterRealPhoton("rhode", Rhode, &PhotoRhode::create);
    RegisterRealPhoton("kokoulin", Kokoulin, &PhotoKokoulin::create);

    RegisterQ2("abramowiczlevinlevymaor91",
             AbramowiczLevinLevyMaor91,
             std::make_pair(&PhotoAbramowiczLevinLevyMaor91::create,
                            &PhotoQ2Interpolant<PhotoAbramowiczLevinLevyMaor91>::create));
    RegisterQ2("abramowiczlevinlevymaor97",
             AbramowiczLevinLevyMaor97,
             std::make_pair(&PhotoAbramowiczLevinLevyMaor97::create,
                            &PhotoQ2Interpolant<PhotoAbramowiczLevinLevyMaor97>::create));
    RegisterQ2("butkevichmikhailov",
             ButkevichMikhailov,
             std::make_pair(&PhotoButkevichMikhailov::create,
                            &PhotoQ2Interpolant<PhotoButkevichMikhailov>::create));
    RegisterQ2("renosarcevicsu",
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

    map_string_to_enum.clear();
}

void PhotonuclearFactory::RegisterShadowEffect(const std::string& name, const Shadow& shadow, RegisterShadowEffectFunction create)
{
    photo_shadow_map_str_[name] = create;
    photo_shadow_map_enum_[shadow] = create;
}

void PhotonuclearFactory::RegisterRealPhoton(const std::string& name, Enum enum_t, RegisterRealPhotonFunction create)
{
    photo_real_map_str_[name] = create;
    photo_real_map_enum_[enum_t] = create;
    map_string_to_enum[name] = enum_t;
}

void PhotonuclearFactory::RegisterQ2(const std::string& name, Enum enum_t, std::pair<RegisterQ2Function, RegisterQ2FunctionInterpolant> create)
{
    photo_q2_map_str_[name] = create;
    photo_q2_map_enum_[enum_t] = create;
    map_string_to_enum[name] = enum_t;
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
// RealPhoton
// ------------------------------------------------------------------------- //

CrossSection* PhotonuclearFactory::CreatePhotoRealPhotonIntegral(const std::string& name,
                                    const ParticleDef& particle_def,
                                    const Medium& medium,
                                    const EnergyCutSettings& cuts,
                                    double multiplier,
                                    bool hardbb) const
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    PhotoRealPhotonMapString::const_iterator it = photo_real_map_str_.find(name_lower);

    if (it != photo_real_map_str_.end())
    {
        RealPhoton* real_photon;

        if (hardbb)
        {
            real_photon = new HardBB(particle_def);
        }
        else
        {
            real_photon = new SoftBB();
        }

        PhotoIntegral* photo = new PhotoIntegral(*it->second(particle_def, medium, cuts, *real_photon, multiplier));
        delete real_photon;
        return photo;
    }
    else
    {
        log_fatal("Photonuclear %s not registerd!", name.c_str());
    }
}

CrossSection* PhotonuclearFactory::CreatePhotoRealPhotonIntegral(const Enum enum_t,
                                    const ParticleDef& particle_def,
                                    const Medium& medium,
                                    const EnergyCutSettings& cuts,
                                    double multiplier,
                                    bool hardbb) const
{
    PhotoRealPhotonMapEnum::const_iterator it = photo_real_map_enum_.find(enum_t);

    if (it != photo_real_map_enum_.end())
    {
        RealPhoton* real_photon;

        if (hardbb)
        {
            real_photon = new HardBB(particle_def);
        }
        else
        {
            real_photon = new SoftBB();
        }

        PhotoIntegral* photo = new PhotoIntegral(*it->second(particle_def, medium, cuts, *real_photon, multiplier));
        delete real_photon;
        return photo;
    }
    else
    {
        log_fatal("Photonuclear %s not registerd!", typeid(enum_t).name());
    }
}

CrossSection* PhotonuclearFactory::CreatePhotoRealPhotonInterpolant(const std::string& name,
                                    const ParticleDef& particle_def,
                                    const Medium& medium,
                                    const EnergyCutSettings& cuts,
                                    double multiplier,
                                    bool hardbb,
                                    InterpolationDef def) const
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    PhotoRealPhotonMapString::const_iterator it = photo_real_map_str_.find(name_lower);

    if (it != photo_real_map_str_.end())
    {
        RealPhoton* real_photon;

        if (hardbb)
        {
            real_photon = new HardBB(particle_def);
        }
        else
        {
            real_photon = new SoftBB();
        }

        PhotoInterpolant* photo = new PhotoInterpolant(*it->second(particle_def, medium, cuts, *real_photon, multiplier), def);
        delete real_photon;
        return photo;
    }
    else
    {
        log_fatal("Photonuclear %s not registerd!", name.c_str());
    }
}

CrossSection* PhotonuclearFactory::CreatePhotoRealPhotonInterpolant(const Enum enum_t,
                                    const ParticleDef& particle_def,
                                    const Medium& medium,
                                    const EnergyCutSettings& cuts,
                                    double multiplier,
                                    bool hardbb,
                                    InterpolationDef def) const
{
    PhotoRealPhotonMapEnum::const_iterator it = photo_real_map_enum_.find(enum_t);

    if (it != photo_real_map_enum_.end())
    {
        RealPhoton* real_photon;

        if (hardbb)
        {
            real_photon = new HardBB(particle_def);
        }
        else
        {
            real_photon = new SoftBB();
        }

        PhotoInterpolant* photo = new PhotoInterpolant(*it->second(particle_def, medium, cuts, *real_photon, multiplier), def);
        delete real_photon;
        return photo;
    }
    else
    {
        log_fatal("Photonuclear %s not registerd!", typeid(enum_t).name());
    }
}

// ------------------------------------------------------------------------- //
// PhotoQ2
// ------------------------------------------------------------------------- //

CrossSection* PhotonuclearFactory::CreatePhotoQ2Integral(const std::string& name,
                                                         const ParticleDef& particle_def,
                                                         const Medium& medium,
                                                         const EnergyCutSettings& cuts,
                                                         double multiplier,
                                                         const Shadow& shadow_effect) const
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    PhotoQ2MapString::const_iterator it = photo_q2_map_str_.find(name_lower);

    if (it != photo_q2_map_str_.end())
    {
        ShadowEffect* shadow = Get().CreateShadowEffect(shadow_effect);

        PhotoIntegral* photo =  new PhotoIntegral(*it->second.first(particle_def, medium, cuts, *shadow, multiplier));
        delete shadow;
        return photo;
    }
    else
    {
        log_fatal("Photonuclear %s not registerd!", name.c_str());
    }
}

CrossSection* PhotonuclearFactory::CreatePhotoQ2Integral(const Enum enum_t,
                                                         const ParticleDef& particle_def,
                                                         const Medium& medium,
                                                         const EnergyCutSettings& cuts,
                                                         double multiplier,
                                                         const Shadow& shadow_effect) const
{
    PhotoQ2MapEnum::const_iterator it = photo_q2_map_enum_.find(enum_t);

    if (it != photo_q2_map_enum_.end())
    {
        ShadowEffect* shadow = Get().CreateShadowEffect(shadow_effect);

        PhotoIntegral* photo =  new PhotoIntegral(*it->second.first(particle_def, medium, cuts, *shadow, multiplier));
        delete shadow;
        return photo;
    }
    else
    {
        log_fatal("Photonuclear %s not registerd!", typeid(enum_t).name());
    }
}

CrossSection* PhotonuclearFactory::CreatePhotoQ2Interpolant(const std::string& name,
                                                            const ParticleDef& particle_def,
                                                            const Medium& medium,
                                                            const EnergyCutSettings& cuts,
                                                            double multiplier,
                                                            const Shadow& shadow_effect,
                                                            InterpolationDef def) const
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    PhotoQ2MapString::const_iterator it = photo_q2_map_str_.find(name_lower);

    if (it != photo_q2_map_str_.end())
    {
        ShadowEffect* shadow = Get().CreateShadowEffect(shadow_effect);

        PhotoInterpolant* photo =  new PhotoInterpolant(*it->second.second(particle_def, medium, cuts, *shadow, multiplier, def), def);
        delete shadow;
        return photo;
    }
    else
    {
        log_fatal("Photonuclear %s not registerd!", name.c_str());
    }
}

CrossSection* PhotonuclearFactory::CreatePhotoQ2Interpolant(const Enum enum_t,
                                                            const ParticleDef& particle_def,
                                                            const Medium& medium,
                                                            const EnergyCutSettings& cuts,
                                                            double multiplier,
                                                            const Shadow& shadow_effect,
                                                            InterpolationDef def) const
{
    PhotoQ2MapEnum::const_iterator it = photo_q2_map_enum_.find(enum_t);

    if (it != photo_q2_map_enum_.end())
    {
        ShadowEffect* shadow = Get().CreateShadowEffect(shadow_effect);

        PhotoInterpolant* photo =  new PhotoInterpolant(*it->second.second(particle_def, medium, cuts, *shadow, multiplier, def), def);
        delete shadow;
        return photo;
    }
    else
    {
        log_fatal("Photonuclear %s not registerd!", typeid(enum_t).name());
    }
}

// ------------------------------------------------------------------------- //
// Most general creator
// ------------------------------------------------------------------------- //

CrossSection* PhotonuclearFactory::CreatePhotonuclear(const Enum enum_t,
                                                      const ParticleDef& particle_def,
                                                      const Medium& medium,
                                                      const EnergyCutSettings& cuts,
                                                      double multiplier,
                                                      const Shadow& shadow_effect,
                                                      bool hardbb,
                                                      bool interpolate,
                                                      InterpolationDef def) const
{
    PhotoQ2MapEnum::const_iterator it_q2 = photo_q2_map_enum_.find(enum_t);
    PhotoRealPhotonMapEnum::const_iterator it_real = photo_real_map_enum_.find(enum_t);

    if (it_q2 != photo_q2_map_enum_.end())
    {
        ShadowEffect* shadow = Get().CreateShadowEffect(shadow_effect);

        if (interpolate)
        {
            PhotoInterpolant* photo =  new PhotoInterpolant(*it_q2->second.second(particle_def, medium, cuts, *shadow, multiplier, def), def);
            delete shadow;
            return photo;
        }
        else
        {
            PhotoIntegral* photo =  new PhotoIntegral(*it_q2->second.first(particle_def, medium, cuts, *shadow, multiplier));
            delete shadow;
            return photo;
        }
    }
    else if (it_real != photo_real_map_enum_.end())
    {
        RealPhoton* real_photon;

        if (hardbb)
        {
            real_photon = new HardBB(particle_def);
        }
        else
        {
            real_photon = new SoftBB();
        }

        if (interpolate)
        {
            PhotoInterpolant* photo = new PhotoInterpolant(*it_real->second(particle_def, medium, cuts, *real_photon, multiplier), def);
            delete real_photon;
            return photo;
        }
        else
        {
            PhotoIntegral* photo = new PhotoIntegral(*it_real->second(particle_def, medium, cuts, *real_photon, multiplier));
            delete real_photon;
            return photo;
        }
    }
    else
    {
        log_fatal("Photonuclear %s not registerd!", typeid(enum_t).name());
    }
}

CrossSection* PhotonuclearFactory::CreatePhotonuclear(const ParticleDef& particle_def,
                                                      const Medium& medium,
                                                      const EnergyCutSettings& cuts,
                                                      const Definition& def,
                                                      bool interpolate,
                                                      InterpolationDef interpolation_def) const
{
    return CreatePhotonuclear(def.parametrization,
                              particle_def,
                              medium,
                              cuts,
                              def.multiplier,
                              def.shadow,
                              def.hardbb,
                              interpolate,
                              interpolation_def);
}

// ------------------------------------------------------------------------- //
PhotonuclearFactory::Enum PhotonuclearFactory::GetEnumFromString(const std::string& name)
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    MapStringToEnum::iterator it = map_string_to_enum.find(name_lower);

    if (it != map_string_to_enum.end())
    {
        return it->second;
    } else
    {
        log_fatal("Photonuclear %s not registerd!", name.c_str());
    }
}
