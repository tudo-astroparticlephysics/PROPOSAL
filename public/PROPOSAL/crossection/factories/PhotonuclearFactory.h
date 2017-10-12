
#pragma once

#include <boost/function.hpp>

#include <map>
#include <string>

#include "PROPOSAL/methods.h"

namespace  PROPOSAL
{

class CrossSection;
class RealPhoton;
class ShadowEffect;
class Photonuclear;
struct ParticleDef;
class EnergyCutSettings;
class Medium;

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

    struct Definition
    {
        Definition()
            : parametrization(AbramowiczLevinLevyMaor97)
            , shadow(ShadowButkevichMikhailov)
            , hardbb(true)
            , multiplier(1.0)
        {
        }

        Enum parametrization;
        Shadow shadow;
        bool hardbb;
        double multiplier;
    };

    // --------------------------------------------------------------------- //
    // Typedefs for readablitiy
    // --------------------------------------------------------------------- //

    typedef boost::function<ShadowEffect*(void)> RegisterShadowEffectFunction;

    typedef boost::function<
        Photonuclear*(const ParticleDef&, const Medium&, const EnergyCutSettings&, const RealPhoton&, double multiplier)>
        RegisterRealPhotonFunction;

    typedef boost::function<
        Photonuclear*(const ParticleDef&, const Medium&, const EnergyCutSettings&, const ShadowEffect&, double multiplier)>
        RegisterQ2Function;

    typedef boost::function<
        Photonuclear*(const ParticleDef&, const Medium&, const EnergyCutSettings&, const ShadowEffect&, double multiplier, InterpolationDef)>
        RegisterQ2FunctionInterpolant;

    typedef std::map<std::string, RegisterShadowEffectFunction > PhotoShadowEffectMapString;
    typedef std::map<Shadow, RegisterShadowEffectFunction > PhotoShadowEffectMapEnum;

    typedef std::map<std::string, RegisterRealPhotonFunction > PhotoRealPhotonMapString;
    typedef std::map<Enum, RegisterRealPhotonFunction > PhotoRealPhotonMapEnum;

    typedef std::map<std::string, std::pair<RegisterQ2Function, RegisterQ2FunctionInterpolant> > PhotoQ2MapString;
    typedef std::map<Enum, std::pair<RegisterQ2Function, RegisterQ2FunctionInterpolant> > PhotoQ2MapEnum;

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

    CrossSection* CreatePhotoRealPhotonIntegral(const std::string&,
                                        const ParticleDef& particle_def,
                                        const Medium& medium,
                                        const EnergyCutSettings& cuts,
                                        double multiplier,
                                        bool hardbb) const;

    CrossSection* CreatePhotoRealPhotonIntegral(const Enum,
                                        const ParticleDef& particle_def,
                                        const Medium& medium,
                                        const EnergyCutSettings& cuts,
                                        double multiplier,
                                        bool hardbb) const;

    CrossSection* CreatePhotoRealPhotonInterpolant(const std::string&,
                                        const ParticleDef& particle_def,
                                        const Medium& medium,
                                        const EnergyCutSettings& cuts,
                                        double multiplier,
                                        bool hardbb,
                                        InterpolationDef = InterpolationDef()) const;

    CrossSection* CreatePhotoRealPhotonInterpolant(const Enum,
                                        const ParticleDef& particle_def,
                                        const Medium& medium,
                                        const EnergyCutSettings& cuts,
                                        double multiplier,
                                        bool hardbb,
                                        InterpolationDef = InterpolationDef()) const;

    // --------------------------------------------------------------------- //
    // Q2 Integration
    // --------------------------------------------------------------------- //

    CrossSection* CreatePhotoQ2Integral(const std::string&,
                                        const ParticleDef&,
                                        const Medium&,
                                        const EnergyCutSettings&,
                                        double multiplier,
                                        const Shadow&) const;

    CrossSection* CreatePhotoQ2Integral(const Enum,
                                        const ParticleDef&,
                                        const Medium&,
                                        const EnergyCutSettings&,
                                        double multiplier,
                                        const Shadow&) const;

    CrossSection* CreatePhotoQ2Interpolant(const std::string&,
                                        const ParticleDef&,
                                        const Medium&,
                                        const EnergyCutSettings&,
                                        double multiplier,
                                        const Shadow&,
                                        InterpolationDef = InterpolationDef()) const;

    CrossSection* CreatePhotoQ2Interpolant(const Enum,
                                        const ParticleDef&,
                                        const Medium&,
                                        const EnergyCutSettings&,
                                        double multiplier,
                                        const Shadow&,
                                        InterpolationDef = InterpolationDef()) const;

    // --------------------------------------------------------------------- //
    // Most general creation
    // --------------------------------------------------------------------- //

    CrossSection* CreatePhotonuclear(const Enum,
                                     const ParticleDef&,
                                     const Medium&,
                                     const EnergyCutSettings&,
                                     double multiplier,
                                     const Shadow&,
                                     bool hardbb,
                                     bool interpolate,
                                     InterpolationDef = InterpolationDef()) const;

    CrossSection* CreatePhotonuclear(const ParticleDef&,
                                     const Medium&,
                                     const EnergyCutSettings&,
                                     const Definition&,
                                     bool interpolate,
                                     InterpolationDef = InterpolationDef()) const;

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
    void RegisterQ2(const std::string& name, Enum, std::pair<RegisterQ2Function, RegisterQ2FunctionInterpolant>);

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

