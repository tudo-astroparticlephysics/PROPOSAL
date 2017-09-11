
#pragma once

#include <boost/function.hpp>

#include <vector>
#include <map>
#include <string>

#include "PROPOSAL/crossection/parametrization/Parametrization.h"

namespace  PROPOSAL
{

class CrossSection;
class RealPhoton;
class ShadowEffect;

class PhotonuclearFactory
{
    public:

    // --------------------------------------------------------------------- //
    // Declare usable enums
    // --------------------------------------------------------------------- //

    enum Enum
    {
        Zeus = 0,
        BezrukovBugaev,
        Rhode,
        Kokoulin,
        AbramowiczLevinLevyMaor91,
        AbramowiczLevinLevyMaor97,
        ButkevichMikhailov,
        RenoSarcevicSu
    };

    enum Shadow
    {
        ShadowDutta,
        ShadowButkevichMikhailov
    };

    // --------------------------------------------------------------------- //
    // Typedefs for readablitiy
    // --------------------------------------------------------------------- //

    typedef boost::function<ShadowEffect*(void)> RegisterShadowEffectFunction;

    typedef boost::function<
        Parametrization*(const ParticleDef&, const Medium&, const EnergyCutSettings&, const RealPhoton&, Parametrization::Definition)>
        RegisterRealPhotonFunction;

    typedef boost::function<
        Parametrization*(const ParticleDef&, const Medium&, const EnergyCutSettings&, const ShadowEffect&, Parametrization::Definition)>
        RegisterQ2Function;

    typedef std::map<std::string, RegisterShadowEffectFunction > PhotoShadowEffectMapString;
    typedef std::map<Shadow, RegisterShadowEffectFunction > PhotoShadowEffectMapEnum;

    typedef std::map<std::string, RegisterRealPhotonFunction > PhotoRealPhotonMapString;
    typedef std::map<Enum, RegisterRealPhotonFunction > PhotoRealPhotonMapEnum;

    typedef std::map<std::string, std::pair<RegisterQ2Function, RegisterQ2Function> > PhotoQ2MapString;
    typedef std::map<Enum, std::pair<RegisterQ2Function, RegisterQ2Function> > PhotoQ2MapEnum;

    typedef std::map<std::string, Enum> MapStringToEnum;

    // --------------------------------------------------------------------- //
    // Create functions
    // --------------------------------------------------------------------- //

    // --------------------------------------------------------------------- //
    // Shadow effect
    // --------------------------------------------------------------------- //

    ShadowEffect* CreateShadowEffect(const std::string&);
    ShadowEffect* CreateShadowEffect(const Shadow&);

    // --------------------------------------------------------------------- //
    // Real Photon
    // --------------------------------------------------------------------- //

    Parametrization* CreatePhotoRealPhotonParam(const std::string&,
                                                const ParticleDef& particle_def,
                                                const Medium& medium,
                                                const EnergyCutSettings& cuts,
                                                bool,
                                                Parametrization::Definition def) const;

    Parametrization* CreatePhotoRealPhotonParam(const Enum,
                                                const ParticleDef& particle_def,
                                                const Medium& medium,
                                                const EnergyCutSettings& cuts,
                                                bool,
                                                Parametrization::Definition def) const;

    CrossSection* CreatePhotoRealPhoton(const std::string&,
                                        const ParticleDef& particle_def,
                                        const Medium& medium,
                                        const EnergyCutSettings& cuts,
                                        bool,
                                        Parametrization::Definition def,
                                        bool interpolate = true) const;

    CrossSection* CreatePhotoRealPhoton(const Enum,
                                        const ParticleDef& particle_def,
                                        const Medium& medium,
                                        const EnergyCutSettings& cuts,
                                        bool,
                                        Parametrization::Definition def,
                                        bool interpolate = true) const;

    // --------------------------------------------------------------------- //
    // Q2 Integration
    // --------------------------------------------------------------------- //

    Parametrization* CreatePhotoQ2IntegralParam(const std::string&,
                                                const ParticleDef&,
                                                const Medium&,
                                                const EnergyCutSettings&,
                                                const Shadow&,
                                                Parametrization::Definition,
                                                bool interpolate = true) const;

    Parametrization* CreatePhotoQ2IntegralParam(const Enum,
                                                const ParticleDef&,
                                                const Medium&,
                                                const EnergyCutSettings&,
                                                const Shadow&,
                                                Parametrization::Definition,
                                                bool interpolate = true) const;

    CrossSection* CreatePhotoQ2Integral(const std::string&,
                                        const ParticleDef&,
                                        const Medium&,
                                        const EnergyCutSettings&,
                                        const Shadow&,
                                        Parametrization::Definition,
                                        bool interpolate = true) const;

    CrossSection* CreatePhotoQ2Integral(const Enum,
                                        const ParticleDef&,
                                        const Medium&,
                                        const EnergyCutSettings&,
                                        const Shadow&,
                                        Parametrization::Definition,
                                        bool interpolate = true) const;

    // --------------------------------------------------------------------- //
    // Most general creation
    // --------------------------------------------------------------------- //

    CrossSection* CreatePhotonuclear(const Enum,
                                     const ParticleDef&,
                                     const Medium&,
                                     const EnergyCutSettings&,
                                     const Shadow&,
                                     bool hardbb,
                                     Parametrization::Definition,
                                     bool interpolate = true) const;

    // ----------------------------------------------------------------------------
    /// @brief Enum string conversation
    // ----------------------------------------------------------------------------
    Enum GetEnumFromString(const std::string&);


    // ----------------------------------------------------------------------------
    /// @brief Register ShadowEffect used for Q2 integration parametrizations
    ///
    /// @param name
    /// @param Enum
    /// @param std::pair
    // ----------------------------------------------------------------------------
    void RegisterShadowEffect(const std::string& name, const Shadow&, RegisterShadowEffectFunction);

    // ----------------------------------------------------------------------------
    /// @brief Register Photonuclear parametrizations
    ///
    /// @param name
    /// @param Enum
    /// @param RegisterRealPhotonFunction
    // ----------------------------------------------------------------------------
    void RegisterRealPhoton(const std::string& name, Enum, RegisterRealPhotonFunction);

    // ----------------------------------------------------------------------------
    /// @brief Register Photonuclear parametrizations
    ///
    /// @param name
    /// @param Enum
    /// @param RegisterQ2Function
    // ----------------------------------------------------------------------------
    void RegisterQ2(const std::string& name, Enum, std::pair<RegisterQ2Function, RegisterQ2Function>);

    // --------------------------------------------------------------------- //
    // Singleton pattern
    // --------------------------------------------------------------------- //

    static PhotonuclearFactory& Get()
    {
        static PhotonuclearFactory instance;
        return instance;
    }

    private:
    PhotonuclearFactory();
    ~PhotonuclearFactory();

    PhotoShadowEffectMapString photo_shadow_map_str_;
    PhotoShadowEffectMapEnum photo_shadow_map_enum_;

    PhotoRealPhotonMapString photo_real_map_str_;
    PhotoRealPhotonMapEnum photo_real_map_enum_;

    PhotoQ2MapString photo_q2_map_str_;
    PhotoQ2MapEnum photo_q2_map_enum_;

    MapStringToEnum map_string_to_enum;
};

} /*  PROPOSAL */

