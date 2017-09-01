/*
 * CrossSections.h
 *
 *  Created on: 2013.03.12
 *      Author: Jan-Hendrik Köhne
 */
#pragma once

// #include <string>

#include "PROPOSAL/particle/PROPOSALParticle.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/EnergyCutSettings.h"

namespace PROPOSAL
{
    class CrossSections;
}

std::ostream& operator<<(std::ostream& os, PROPOSAL::CrossSections const&crossSections);

namespace PROPOSAL{

namespace ParametrizationType
{
    enum Enum
    {
        BremsKelnerKokoulinPetrukhin  = 11,
        BremsAndreevBezrukovBugaev    = 12,
        BremsPetrukhinShestakov       = 13,
        BremsCompleteScreeningCase    = 14,

        PhotoKokoulinShadowBezrukovSoft                = 31,
        PhotoKokoulinShadowBezrukovHard                = -31,
        PhotoRhodeShadowBezrukovSoft                   = 32,
        PhotoRhodeShadowBezrukovHard                   = -32,
        PhotoBezrukovBugaevShadowBezrukovSoft          = 33,
        PhotoBezrukovBugaevShadowBezrukovHard          = -33,
        PhotoZeusShadowBezrukovSoft                    = 34,
        PhotoZeusShadowBezrukovHard                    = -34,
        PhotoAbramowiczLevinLevyMaor91ShadowDutta      = 35,
        PhotoAbramowiczLevinLevyMaor91ShadowButkevich  = -35,
        PhotoAbramowiczLevinLevyMaor97ShadowDutta      = 36,
        PhotoAbramowiczLevinLevyMaor97ShadowButkevich  = -36,
        PhotoButkevichMikhailovShadowDutta             = 37,
        PhotoButkevichMikhailovShadowButkevich         = -37,
        PhotoRenoSarcevicSuShadowDutta                 = 38,
        PhotoRenoSarcevicSuShadowButkevich             = -38,

        EPairKelnerKokoulinPetrukhin  = 51,

        IonizBetheBloch  = 71
    };
}

/*! \class CrossSections CrossSections.h "CrossSections.h"
 *  \brief This is a pure virtual class
 */
class CrossSections
{


protected:

    std::string name_;
    ParticleType::Enum type_;

    // PROPOSALParticle*   particle_;
    // PROPOSALParticle*   backup_particle_;
    Medium*     medium_;
    EnergyCutSettings* cut_settings_;


    double ebig_;   //!< upper bound of parameterizations

    double rnd_;    //!< This random number will be stored in CalculateDNdx to avoid calculate dNdx a second time in ClaculateSochasticLoss when it is already done


    // Interpolation flags
    bool do_dedx_Interpolation_;
    bool do_dndx_Interpolation_;

    //CrossSection multiplier
    double multiplier_;

    ParametrizationType::Enum  parametrization_;

    bool    lpm_effect_enabled_;
    bool    init_lpm_effect_;

    int     order_of_interpolation_;
    double  sum_of_rates_;

//----------------------------------------------------------------------------//

    virtual double FunctionToDEdxIntegral(const PROPOSALParticle&, double variable) = 0;

//----------------------------------------------------------------------------//

    virtual double CalculateStochasticLoss(const PROPOSALParticle&, double rnd1) = 0;

//----------------------------------------------------------------------------//

public:

    //bounds of integration
    struct IntegralLimits
    {
        double vMax; //!< upper bound of integration
        double vUp;  //!< lower bound of integration
        double vMin; //!< lowest physical possible bound of integration
    };

    //Constructor
    CrossSections();
//----------------------------------------------------------------------------//

    CrossSections(Medium* medium, EnergyCutSettings* cut_settings);
    CrossSections(const CrossSections& crossSections);
    bool operator==(const CrossSections &crossSections) const;
    bool operator!=(const CrossSections &crossSections) const;
    friend std::ostream& operator<<(std::ostream& os, CrossSections const&crossSections);

//----------------------------------------------------------------------------//

    // Memberfunctions
    void swap(CrossSections &crossSections);


//----------------------------------------------------------------------------//

    virtual double CalculatedEdx(const PROPOSALParticle&) = 0;

//----------------------------------------------------------------------------//

    virtual double CalculatedNdx(const PROPOSALParticle&) = 0;

//----------------------------------------------------------------------------//

    virtual double CalculatedNdx(const PROPOSALParticle&, double rnd) = 0;

//----------------------------------------------------------------------------//

    virtual double CalculateStochasticLoss(const PROPOSALParticle&, double rnd1, double rnd2) = 0;

//----------------------------------------------------------------------------//

    virtual void EnableDNdxInterpolation(const PROPOSALParticle&, std::string path ="", bool raw=false) = 0;

//----------------------------------------------------------------------------//

    virtual void EnableDEdxInterpolation(const PROPOSALParticle&, std::string path ="", bool raw=false) = 0;

//----------------------------------------------------------------------------//

    virtual void DisableDNdxInterpolation() = 0;

//----------------------------------------------------------------------------//

    virtual void DisableDEdxInterpolation() = 0;

//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//

    virtual double FunctionToDNdxIntegral(const PROPOSALParticle&, double variable) = 0;

//----------------------------------------------------------------------------//

    virtual CrossSections::IntegralLimits SetIntegralLimits(const PROPOSALParticle&, int component) = 0;


//----------------------------------------------------------------------------//
//-----------------------------Getter and Setter------------------------------//
//----------------------------------------------------------------------------//
    // Setter

    void SetParametrizationLimit(const double ebig){ ebig_ = ebig; }

    void SetMultiplier(const double multiplier){ multiplier_ = multiplier; }

    void SetMedium(Medium* medium){ medium_ = medium; }

    void SetEnergyCutSettings(EnergyCutSettings* cuts){ cut_settings_ = cuts; }

    void EnableLpmEffect(const bool lpm_effect_enabled){ lpm_effect_enabled_ = lpm_effect_enabled; }


    // virtual void ValidateOptions() = 0;

//----------------------------------------------------------------------------//
    // Getter

    std::string GetName() const { return name_; }

    ParticleType::Enum GetType() const { return type_; }

    double GetEbig() const { return ebig_; }

    Medium* GetMedium() const { return medium_; }

    double GetMultiplier() const { return multiplier_; }

    ParametrizationType::Enum GetParametrization() const { return parametrization_; }

    bool GetLpmEffectEnabled() const { return lpm_effect_enabled_; }

    EnergyCutSettings* GetEnergyCutSettings() const { return cut_settings_; }

//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//


    // destructor

    virtual ~CrossSections(){}

    PROPOSALParticle *GetBackup_particle() const;
    void SetBackup_particle(PROPOSALParticle *backup_particle);
    void RestoreBackup_particle();
};

}
