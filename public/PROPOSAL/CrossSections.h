/*
 * CrossSections.h
 *
 *  Created on: 2013.03.12
 *      Author: Jan-Hendrik Köhne
 */

#ifndef CrossSections_H
#define CrossSections_H


#include "PROPOSAL/MathModel.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Particle.h"
#include "PROPOSAL/Medium.h"
#include "PROPOSAL/EnergyCutSettings.h"
#include <string>
#include <boost/program_options.hpp>



/*! \class CrossSections CrossSections.h "CrossSections.h"
    \brief This is a pure virtual class
 */



class CrossSections
{


protected:

    std::string name_;

    Particle*   particle_;
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

    int         parametrization_;

    bool        lpm_effect_enabled_;
    bool        init_lpm_effect_;

    int         order_of_interpolation_;
    double      sum_of_rates_;
//----------------------------------------------------------------------------//

    virtual double FunctionToDEdxIntegral(double variable) = 0;

//----------------------------------------------------------------------------//

    virtual double CalculateStochasticLoss(double rnd1) = 0;

//----------------------------------------------------------------------------//

public:

    //Constructor
    CrossSections();
//----------------------------------------------------------------------------//

    CrossSections(Particle* particle, Medium* medium, EnergyCutSettings* cut_settings);
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

    void SetParticle(Particle *particle);

//----------------------------------------------------------------------------//

    void SetMedium(Medium *medium);

//----------------------------------------------------------------------------//

    void SetVMin(double vMin=0);

//----------------------------------------------------------------------------//

    void SetVMax(double vMax=0);

//----------------------------------------------------------------------------//

    void SetVUp(double vUp=0);

//----------------------------------------------------------------------------//

    virtual void SetParametrization(int parametrization=1) = 0;

//----------------------------------------------------------------------------//

    void SetEnergyCutSettings(EnergyCutSettings *cuts);

//----------------------------------------------------------------------------//

    virtual boost::program_options::options_description CreateOptions() = 0;

//----------------------------------------------------------------------------//
    virtual void ValidateOptions() = 0;

//----------------------------------------------------------------------------//
    // Getter
    std::string GetName() const
    {
        return name_;
    }
//----------------------------------------------------------------------------//
    double GetEbig() const
    {
        return ebig_;
    }

//----------------------------------------------------------------------------//

    Particle* GetParticle() const
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

    double GetParametrization() const
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

    void SetShadow(int newshadow);
//----------------------------------------------------------------------------//

    // destructor

    virtual ~CrossSections(){}

};



#endif //CrossSections_H
