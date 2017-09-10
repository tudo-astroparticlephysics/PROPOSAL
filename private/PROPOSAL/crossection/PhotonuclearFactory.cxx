
#include <boost/algorithm/string.hpp> // case insensitive string compare for configuration file

#include "PROPOSAL/crossection/PhotoIntegral.h"
#include "PROPOSAL/crossection/PhotoInterpolant.h"
#include "PROPOSAL/crossection/parametrization/Photonuclear.h"
#include "PROPOSAL/crossection/parametrization/PhotoRealPhotonAssumption.h"
#include "PROPOSAL/crossection/parametrization/PhotoQ2Integration.h"
#include "PROPOSAL/crossection/PhotonuclearFactory.h"

#include "PROPOSAL/Output.h"

using namespace PROPOSAL;

PhotonuclearFactory::PhotonuclearFactory()
{
    // Register all photonuclear parametrizations in lower case!

    Register("Zeus", Zeus, &PhotoZeus::create);
    Register("BezrukovBugaev", BezrukovBugaev, &PhotoBezrukovBugaev::create);
    Register("Rhode", Rhode, &PhotoRhode::create);
    Register("Kokoulin", Kokoulin, &PhotoKokoulin::create);
    RegisterQ2("AbramowiczLevinLevyMaor91",
             AbramowiczLevinLevyMaor91,
             std::make_pair(&PhotoAbramowiczLevinLevyMaor91::create,
                            &PhotoQ2Interpolant<PhotoAbramowiczLevinLevyMaor91>::create));
    RegisterQ2("AbramowiczLevinLevyMaor97",
             AbramowiczLevinLevyMaor97,
             std::make_pair(&PhotoAbramowiczLevinLevyMaor97::create,
                            &PhotoQ2Interpolant<PhotoAbramowiczLevinLevyMaor97>::create));
    RegisterQ2("ButkevichMikhailov",
             ButkevichMikhailov,
             std::make_pair(&PhotoButkevichMikhailov::create,
                            &PhotoQ2Interpolant<PhotoButkevichMikhailov>::create));
    RegisterQ2("RenoSarcevicSu",
             RenoSarcevicSu,
             std::make_pair(&PhotoRenoSarcevicSu::create,
                            &PhotoQ2Interpolant<PhotoRenoSarcevicSu>::create));
    // Register("AbramowiczLevinLevyMaor97", AbramowiczLevinLevyMaor97, &PhotoAbramowiczLevinLevyMaor97::create);
    // Register("ButkevichMikhailov", ButkevichMikhailov, &PhotoButkevichMikhailov::create);
    // Register("RenoSarcevicSu", RenoSarcevicSu, &PhotoRenoSarcevicSu::create);
}

PhotonuclearFactory::~PhotonuclearFactory()
{
    photo_real_map_str_.clear();
    photo_real_map_enum_.clear();
    photo_q2_map_str_.clear();
    photo_q2_map_enum_.clear();
    map_string_to_enum.clear();
}

void PhotonuclearFactory::Register(const std::string& name, Enum enum_t, RegisterRealPhotonFunction create)
{
    photo_real_map_str_[name] = create;
    photo_real_map_enum_[enum_t] = create;
    map_string_to_enum[name] = enum_t;
}

void PhotonuclearFactory::RegisterQ2(const std::string& name, Enum enum_t, std::pair<RegisterQ2Function, RegisterQ2Function> create)
{
    photo_q2_map_str_[name] = create;
    photo_q2_map_enum_[enum_t] = create;
    map_string_to_enum[name] = enum_t;
}

// ------------------------------------------------------------------------- //
Parametrization* PhotonuclearFactory::CreatePhotoRealPhotonParam(const std::string& name,
                                                                 const ParticleDef& particle_def,
                                                                 const Medium& medium,
                                                                 const EnergyCutSettings& cuts,
                                                                 const RealPhoton& real_photon,
                                                                 Parametrization::Definition def) const
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    PhotoRealPhotonMapString::const_iterator it = photo_real_map_str_.find(name_lower);

    if (it != photo_real_map_str_.end())
    {
        return it->second(particle_def, medium, cuts, real_photon, def);
    } else
    {
        log_fatal("Photonuclear %s not registerd!", name.c_str());
    }
}

// ------------------------------------------------------------------------- //
Parametrization* PhotonuclearFactory::CreatePhotoRealPhotonParam(Enum enum_t,
                                                                 const ParticleDef& particle_def,
                                                                 const Medium& medium,
                                                                 const EnergyCutSettings& cuts,
                                                                 const RealPhoton& real_photon,
                                                                 Parametrization::Definition def) const
{
    PhotoRealPhotonMapEnum::const_iterator it = photo_real_map_enum_.find(enum_t);

    if (it != photo_real_map_enum_.end())
    {
        return it->second(particle_def, medium, cuts, real_photon, def);
    } else
    {
        log_fatal("Photonuclear %s not registerd!", typeid(enum_t).name());
    }
}

// ------------------------------------------------------------------------- //
CrossSection* PhotonuclearFactory::CreatePhotoRealPhoton(const std::string& name,
                                                      const ParticleDef& particle_def,
                                                      const Medium& medium,
                                                      const EnergyCutSettings& cuts,
                                                      const RealPhoton& real_photon,
                                                      Parametrization::Definition def,
                                                      bool interpolate) const
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    PhotoRealPhotonMapString::const_iterator it = photo_real_map_str_.find(name_lower);

    if (it != photo_real_map_str_.end())
    {
        if (interpolate)
        {
            return new PhotoInterpolant(*it->second(particle_def, medium, cuts, real_photon, def));
        }
        else
        {
            return new PhotoIntegral(*it->second(particle_def, medium, cuts, real_photon, def));
        }
    } else
    {
        log_fatal("Photonuclear %s not registerd!", name.c_str());
    }
}

// ------------------------------------------------------------------------- //
CrossSection* PhotonuclearFactory::CreatePhotoRealPhoton(const Enum enum_t,
                                                      const ParticleDef& particle_def,
                                                      const Medium& medium,
                                                      const EnergyCutSettings& cuts,
                                                      const RealPhoton& real_photon,
                                                      Parametrization::Definition def,
                                                      bool interpolate) const
{
    PhotoRealPhotonMapEnum::const_iterator it = photo_real_map_enum_.find(enum_t);

    if (it != photo_real_map_enum_.end())
    {
        if (interpolate)
        {
            return new PhotoInterpolant(*it->second(particle_def, medium, cuts, real_photon, def));
        }
        else
        {
            return new PhotoIntegral(*it->second(particle_def, medium, cuts, real_photon, def));
        }
    } else
    {
        log_fatal("Photonuclear %s not registerd!", typeid(enum_t).name());
    }
}

// ------------------------------------------------------------------------- //
Parametrization* PhotonuclearFactory::CreatePhotoQ2IntegralParam(const std::string& name,
                                                                 const ParticleDef& particle_def,
                                                                 const Medium& medium,
                                                                 const EnergyCutSettings& cuts,
                                                                 const ShadowEffect& shadow_effect,
                                                                 Parametrization::Definition def,
                                                                 bool interpolate) const
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    PhotoQ2MapString::const_iterator it = photo_q2_map_str_.find(name_lower);

    if (it != photo_q2_map_str_.end())
    {
        if (interpolate)
        {
            return it->second.second(particle_def, medium, cuts, shadow_effect, def);
        }
        else
        {
            return it->second.first(particle_def, medium, cuts, shadow_effect, def);
        }
    } else
    {
        log_fatal("Photonuclear %s not registerd!", name.c_str());
    }
}

// ------------------------------------------------------------------------- //
Parametrization* PhotonuclearFactory::CreatePhotoQ2IntegralParam(const Enum enum_t,
                                                                 const ParticleDef& particle_def,
                                                                 const Medium& medium,
                                                                 const EnergyCutSettings& cuts,
                                                                 const ShadowEffect& shadow_effect,
                                                                 Parametrization::Definition def,
                                                                 bool interpolate) const
{
    PhotoQ2MapEnum::const_iterator it = photo_q2_map_enum_.find(enum_t);

    if (it != photo_q2_map_enum_.end())
    {
        if (interpolate)
        {
            return it->second.second(particle_def, medium, cuts, shadow_effect, def);
        }
        else
        {
            return it->second.first(particle_def, medium, cuts, shadow_effect, def);
        }
    } else
    {
        log_fatal("Photonuclear %s not registerd!", typeid(enum_t).name());
    }
}

// ------------------------------------------------------------------------- //
CrossSection* PhotonuclearFactory::CreatePhotoQ2Integral(const std::string& name,
                                                         const ParticleDef& particle_def,
                                                         const Medium& medium,
                                                         const EnergyCutSettings& cuts,
                                                         const ShadowEffect& shadow_effect,
                                                         Parametrization::Definition def,
                                                         bool interpolate) const
{
    std::string name_lower = boost::algorithm::to_lower_copy(name);

    PhotoQ2MapString::const_iterator it = photo_q2_map_str_.find(name_lower);

    if (it != photo_q2_map_str_.end())
    {
        if (interpolate)
        {
            return new PhotoInterpolant(*it->second.second(particle_def, medium, cuts, shadow_effect, def));
        }
        else
        {
            return new PhotoIntegral(*it->second.first(particle_def, medium, cuts, shadow_effect, def));
        }
    } else
    {
        log_fatal("Photonuclear %s not registerd!", name.c_str());
    }
}

// ------------------------------------------------------------------------- //
CrossSection* PhotonuclearFactory::CreatePhotoQ2Integral(const Enum enum_t,
                                                         const ParticleDef& particle_def,
                                                         const Medium& medium,
                                                         const EnergyCutSettings& cuts,
                                                         const ShadowEffect& shadow_effect,
                                                         Parametrization::Definition def,
                                                         bool interpolate) const
{
    PhotoQ2MapEnum::const_iterator it = photo_q2_map_enum_.find(enum_t);

    if (it != photo_q2_map_enum_.end())
    {
        if (interpolate)
        {
            return new PhotoInterpolant(*it->second.second(particle_def, medium, cuts, shadow_effect, def));
        }
        else
        {
            return new PhotoIntegral(*it->second.first(particle_def, medium, cuts, shadow_effect, def));
        }
    } else
    {
        log_fatal("Photonuclear %s not registerd!", typeid(enum_t).name());
    }
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
