
/******************************************************************************
 *																			  *
 * This file is part of the simulation tool PROPOSAL.						  *
 *																			  *
 * Copyright (C) 2017 TU Dortmund University, Department of Physics,		  *
 * 				      Chair Experimental Physics 5b							  *
 *																			  *
 * This software may be modified and distributed under the terms of a		  *
 * modified GNU Lesser General Public Licence version 3 (LGPL),               *
 * copied verbatim in the file "LICENSE". 									  *
 *																			  *
 * Modifcations to the LGPL License:										  *
 *																			  *
 *      1. The user shall acknowledge the use of PROPOSAL by citing the		  *
 *         following reference:												  *
 *																			  *
 *         J.H. Koehne et al.  Comput.Phys.Commun. 184 (2013) 2070-2090 DOI:  *
 *         10.1016/j.cpc.2013.04.001										  *
 *																			  *
 *      2. The user should report any bugs/errors or improvments to the       *
 *         current maintainer of PROPOSAL or open an issue on the             *
 *		   GitHub webpage													  *
 *																			  *
 *		   "https://github.com/tudo-astroparticlephysics/PROPOSAL"			  *
 *																			  *
 ******************************************************************************/

#pragma once

#ifndef CrossSections_H
#define CrossSections_H


// #include <string>

#include "PROPOSAL/MathModel.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/PROPOSALParticle.h"
#include "PROPOSAL/Medium.h"
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

    PROPOSALParticle*   particle_;
    PROPOSALParticle*   backup_particle_;
    Medium*     medium_;
    EnergyCutSettings* cut_settings_;


    //bounds of integration
    double vMax_;   //!< upper bound of integration
    double vUp_;    //!< lower bound of integration
    double vMin_;   //!< lowest physical possible bound of integration

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

    virtual double FunctionToDEdxIntegral(double variable) = 0;

//----------------------------------------------------------------------------//

    virtual double CalculateStochasticLoss(double rnd1) = 0;

//----------------------------------------------------------------------------//

public:

    //Constructor
    CrossSections();
//----------------------------------------------------------------------------//

    CrossSections(PROPOSALParticle* particle, Medium* medium, EnergyCutSettings* cut_settings);
    CrossSections(const CrossSections& crossSections);
    bool operator==(const CrossSections &crossSections) const;
    bool operator!=(const CrossSections &crossSections) const;
    friend std::ostream& operator<<(std::ostream& os, CrossSections const&crossSections);

//----------------------------------------------------------------------------//

    // Memberfunctions
    void swap(CrossSections &crossSections);


//----------------------------------------------------------------------------//

    virtual double CalculatedEdx() = 0;

//----------------------------------------------------------------------------//

    virtual double CalculatedNdx() = 0;

//----------------------------------------------------------------------------//

    virtual double CalculatedNdx(double rnd) = 0;

//----------------------------------------------------------------------------//

    virtual double CalculateStochasticLoss(double rnd1, double rnd2) = 0;

//----------------------------------------------------------------------------//

    virtual void EnableDNdxInterpolation(std::string path ="", bool raw=false) = 0;

//----------------------------------------------------------------------------//

    virtual void EnableDEdxInterpolation(std::string path ="", bool raw=false) = 0;

//----------------------------------------------------------------------------//

    virtual void DisableDNdxInterpolation() = 0;

//----------------------------------------------------------------------------//

    virtual void DisableDEdxInterpolation() = 0;

//----------------------------------------------------------------------------//

    void SetParametrizationLimit(double ebig=BIGENERGY);

//----------------------------------------------------------------------------//

    virtual double FunctionToDNdxIntegral(double variable) = 0;

//----------------------------------------------------------------------------//

    virtual void SetIntegralLimits(int component) = 0;

//----------------------------------------------------------------------------//

    // Setter

//----------------------------------------------------------------------------//
    void SetMultiplier(double multiplier=1.);

//----------------------------------------------------------------------------//

    void SetParticle(PROPOSALParticle *particle);

//----------------------------------------------------------------------------//

    void SetMedium(Medium *medium);

//----------------------------------------------------------------------------//

    void SetVMin(double vMin=0);

//----------------------------------------------------------------------------//

    void SetVMax(double vMax=0);

//----------------------------------------------------------------------------//

    void SetVUp(double vUp=0);

//----------------------------------------------------------------------------//

    virtual void SetParametrization(
        ParametrizationType::Enum parametrization = ParametrizationType::BremsKelnerKokoulinPetrukhin) = 0;

//----------------------------------------------------------------------------//

    void SetEnergyCutSettings(EnergyCutSettings *cuts);

//----------------------------------------------------------------------------//

    virtual void ValidateOptions() = 0;

//----------------------------------------------------------------------------//
    // Getter
    std::string GetName() const
    {
        return name_;
    }

    ParticleType::Enum GetType() const
    {
        return type_;
    }
//----------------------------------------------------------------------------//
    double GetEbig() const
    {
        return ebig_;
    }

//----------------------------------------------------------------------------//

    PROPOSALParticle* GetParticle() const
    {
        return particle_;
    }

//----------------------------------------------------------------------------//
    Medium* GetMedium() const
    {
        return medium_;
    }

//----------------------------------------------------------------------------//
    double GetMultiplier() const
    {
        return multiplier_;
    }

//----------------------------------------------------------------------------//

    double GetVMax() const
    {
        return vMax_;
    }

//----------------------------------------------------------------------------//

    double GetVMin() const
    {
        return vMin_;
    }

//----------------------------------------------------------------------------//

    double GetVUp() const
    {
        return vUp_;
    }

//----------------------------------------------------------------------------//

    ParametrizationType::Enum GetParametrization() const
    {
        return parametrization_;
    }

//----------------------------------------------------------------------------//

    bool GetLpmEffectEnabled() const
    {
        return lpm_effect_enabled_;
    }
//----------------------------------------------------------------------------//
    EnergyCutSettings* GetEnergyCutSettings() const
    {
        return cut_settings_;
    }

//----------------------------------------------------------------------------//

    void EnableLpmEffect(bool lpm_effect_enabled);

//----------------------------------------------------------------------------//


    // destructor

    virtual ~CrossSections(){}

    PROPOSALParticle *GetBackup_particle() const;
    void SetBackup_particle(PROPOSALParticle *backup_particle);
    void RestoreBackup_particle();
};

}


#endif //CrossSections_H
