
/******************************************************************************
 *                                                                            *
 * This file is part of the simulation tool PROPOSAL.                         *
 *                                                                            *
 * Copyright (C) 2017 TU Dortmund University, Department of Physics,          *
 *                    Chair Experimental Physics 5b                           *
 *                                                                            *
 * This software may be modified and distributed under the terms of a         *
 * modified GNU Lesser General Public Licence version 3 (LGPL),               *
 * copied verbatim in the file "LICENSE".                                     *
 *                                                                            *
 * Modifcations to the LGPL License:                                          *
 *                                                                            *
 *      1. The user shall acknowledge the use of PROPOSAL by citing the       *
 *         following reference:                                               *
 *                                                                            *
 *         J.H. Koehne et al.  Comput.Phys.Commun. 184 (2013) 2070-2090 DOI:  *
 *         10.1016/j.cpc.2013.04.001                                          *
 *                                                                            *
 *      2. The user should report any bugs/errors or improvments to the       *
 *         current maintainer of PROPOSAL or open an issue on the             *
 *         GitHub webpage                                                     *
 *                                                                            *
 *         "https://github.com/tudo-astroparticlephysics/PROPOSAL"            *
 *                                                                            *
 ******************************************************************************/

#pragma once

#ifndef ProcessCollection_H
#define ProcessCollection_H

// #include <vector>
// #include <utility>

#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Epairproduction.h"
#include "PROPOSAL/Ionization.h"
#include "PROPOSAL/Photonuclear.h"
#include "PROPOSAL/Decay.h"
#include "PROPOSAL/ContinuousRandomization.h"
#include "PROPOSAL/Geometry.h"
#include "PROPOSAL/Scattering.h"
// #include "PROPOSAL/Constants.h"
// #include "PROPOSAL/EnergyCutSettings.h"
// #include "PROPOSAL/PROPOSALParticle.h"
// #include "PROPOSAL/Integral.h"
// #include "PROPOSAL/Interpolant.h"
// #include "PROPOSAL/Medium.h"
// #include "PROPOSAL/CrossSections.h"

namespace PROPOSAL
{
    class ProcessCollection;
}

std::ostream& operator<<(std::ostream& os, PROPOSAL::ProcessCollection const& collection);

namespace PROPOSAL{

/*! \class ProcessCollection ProcessCollection.h "CrossSections.h"
    \brief initializes all cross sections and keeps references to them
 */

class ProcessCollection : public MathModel
{

protected:

    int         order_of_interpolation_;
    bool        do_interpolation_;
    bool        lpm_effect_enabled_;
    double      ini_;
    bool        debug_;
    bool        do_weighting_;           //!< Do weigthing? Set to false in constructor
    double      weighting_order_;        //!< Re-weighting order. Set to 0 in constructor
    double      weighting_starts_at_;    //!< Distance at which re-weighting starts. Set to 0 in constructor

    bool        enable_randomization_;        //!< if true continuous randomization will be enabled (to remember if randomization should be enable when cross sections are initalized)
    bool        do_continuous_randomization_; //!< if true randomization of continuous energy losses is enabled
    bool        do_scattering_;               //!< if true moliere scattering is enabled
    int         location_;                    //!< 0 = infront of the detector, 1 = inside the detector, 2 = behind the detector

    double      density_correction_;          //!< density correction factor

    bool    do_time_interpolation_;     //!< If true, CalculateParticleTime is interpolated
    bool    do_exact_time_calulation_;  //!< exact local time calculation enabled if true


    Geometry*   geometry_;

    /*!
     * \brief indicates if the interpolated function is increasing or decreasing.
     *
     * True  =   Interpolated function is increasing. \n
     * False =   Interpolated function is decreasing. BigLow is set to f(x_min)
     *           and the interpolation points that are saved are BigLow-f(x).
    */
    bool up_;

    std::vector<double> bigLow_;               //!< See the describtion of the variable up.
    std::vector<double> storeDif_;             //!< Stores the interpolated values for farther calculations

    Interpolant*        interpolant_;
    Interpolant*        interpolant_diff_;
    PROPOSALParticle*           particle_;
    PROPOSALParticle*           backup_particle_;
    Medium*             medium_;
    Integral*           integral_;
    EnergyCutSettings*  cut_settings_;


    std::vector<CrossSections*> crosssections_;
    Decay* decay_;
    Integral* prop_decay_;
    Integral* prop_interaction_;
    Interpolant *interpol_prop_decay_;           //!< Interpolant object of the Integral of the function FunctionToPropIntegralDecay
    Interpolant *interpol_prop_decay_diff_;      //!< Interpolant object of the function FunctionToPropIntegralDecay
    Interpolant *interpol_prop_interaction_;     //!< Interpolant object of the Integral of the function FunctionToPropIntegralInteraction
    Interpolant *interpol_prop_interaction_diff_;//!< Interpolant object of the function FunctionToPropIntegralInteraction

    ContinuousRandomization *randomizer_;
    Scattering* scattering_;

    Interpolant* interpol_time_particle_;
    Interpolant* interpol_time_particle_diff_;

    Integral*    time_particle_;

    //Memeberfunctions

    /**
     * function for energy range calculation - interface to Integral
     * For decay
     *
     *  \param  energy particle energy
     *  \return Returns the probability [1/MeV]
     */

    double FunctionToPropIntegralDecay(double energy);

//----------------------------------------------------------------------------//

    /**
     * function for energy range calculation - interface to Integral
     * For Interaction
     *
     *  \param  energy particle energy
     *  \return Returns the probability [1/MeV]
     */

    double FunctionToPropIntegralInteraction(double energy);

//----------------------------------------------------------------------------//

    /*!
    * function for time delta calculation - interface to Integral
    *
    */

    double FunctionToTimeIntegral(double E);

//----------------------------------------------------------------------------//

    double InterpolTimeParticle(double energy);

//----------------------------------------------------------------------------//
    double InterpolTimeParticleDiff(double energy);

//----------------------------------------------------------------------------//

    /*!
    1d parametrization ;
    \f[return= \int_{e}^{e_{low}} f(E) dE \f] with \f$
    e_{low} \f$ energy below which the particle is lost
    \param e energy [MeV]
    */
    double FunctionToBuildInterpolant(double energy);

//----------------------------------------------------------------------------//

    /*!
    1d parametrization ;
    \f[return=function(e)\f];
    \param e energy [MeV]
    */
    double FunctionToBuildInterpolantDiff(double energy);

//----------------------------------------------------------------------------//

    double InterpolPropDecay(double energy);

//----------------------------------------------------------------------------//
    double InterpolPropDecayDiff(double energy);

//----------------------------------------------------------------------------//
    double InterpolPropInteraction(double energy);

//----------------------------------------------------------------------------//
    double InterpolPropInteractionDiff(double energy);


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
public:

    // constructors

    ProcessCollection();
    ProcessCollection(const ProcessCollection& collection);
    ProcessCollection& operator=(const ProcessCollection& collection);
    bool operator==(const ProcessCollection &collection) const;
    bool operator!=(const ProcessCollection &collection) const;
    friend std::ostream& operator<<(std::ostream& os, ProcessCollection const& collection);

//----------------------------------------------------------------------------//

    /// @brief  initializes all cross sections,

    ProcessCollection(PROPOSALParticle *particle, Medium *medium, EnergyCutSettings* cut_settings);

//----------------------------------------------------------------------------//

    // Memberfunctions
    void swap(ProcessCollection &collection);

//----------------------------------------------------------------------------//
    /*!
    * time delta, corresponding to the given propagation distance
    *
    * \param    ei  initial energy
    * \param    ef  final energy
    * \return   time delta
    */

    double CalculateParticleTime(double ei, double ef);

//----------------------------------------------------------------------------//

    /*!
    returns the value of the distance integral from ei to ef;
    \f[ dx = \int_{e_i}^{e_f} - \frac{1}{  \frac{dE}{dx}\big|_{Ioniz}
    +\frac{dE}{dx}\big|_{Brems}+\frac{dE}{dx}\big|_{PhotoNuc}+\frac{dE}{dx}
    \big|_{epair}  } dE \f]
    \return value of the distance integral from ei to ef
    \param ei lower integration limit (initial energy)
    \param ef upper integration limit (final energy)
    \param dist ???
    */
    double CalculateDisplacement(double ei, double ef, double dist);

//----------------------------------------------------------------------------//

    /*!
    final energy, corresponding to the given value of displacement dist;
    returns \f$e_f\f$ which fullfills
    \f[\int_{e_i}^{\hat{e}_f} - \frac{1}{  \frac{dE}{dx}\big|_{Ioniz} +
    \frac{dE}{dx}\big|_{Brems}+\frac{dE}{dx}\big|_{PhotoNuc}+\frac{dE}{dx}
    \big|_{epair}  } dE = dist\f]
    \return final energy [MeV]
    \param ei initial energy [MeV]
    \param dist value of displacement
    */

    double CalculateFinalEnergy(double ei, double dist);

//----------------------------------------------------------------------------//

    /**
     * final energy, corresponding to the given value rnd of the
     * tracking integral
     *
     *  \param  ei                      initial energy
     *  \param  rnd                     random number which is used for the calculation
     *  \param  particle_interaction    particle interaction? (false = decay)
     *  \return final energy due to continous energy losses [MeV]
     */

    double CalculateFinalEnergy(double ei, double rnd, bool particle_interaction);

//----------------------------------------------------------------------------//
    /**
     * Calculates the value of the tracking integral
     *
     *  \param  initial_energy          initial energy
     *  \param  rnd                     random number which is used for the calculation
     *  \param  particle_interaction    particle interaction? (false = decay)
     *  \return value of the tracking integral [ 1 ]
     */

    double CalculateTrackingIntegal(double initial_energy, double rnd, bool particle_interaction);

//----------------------------------------------------------------------------//
    /**
     * randomize the continuous energy loss
     *
     *  \param  initial_energy          initial energy
     *  \param  final_energy            final energy
     */

    double Randomize(double initial_energy, double final_energy);

//----------------------------------------------------------------------------//

    /**
     *  Makes Stochastic Energyloss
     *
     *  \return pair of energy loss [MeV] and kind of interaction
     */
    std::pair<double, ParticleType::Enum> MakeStochasticLoss();

    std::pair<double, ParticleType::Enum> MakeStochasticLoss(double rnd1,double rnd2, double rnd3);
//----------------------------------------------------------------------------//

    /**
     *  Makes Decay
     *
     *  \return pair of product energy [MeV] and kind of product
     */
    std::pair<double, ParticleType::Enum> MakeDecay();

    std::pair<double, ParticleType::Enum> MakeDecay(double rnd1,double rnd2, double rnd3);
//----------------------------------------------------------------------------//
    /**
     * Enables the Interpolation including dEdx and dNdx for
     * every crosssection in vector crosssections_
     */
    void EnableInterpolation(std::string path ="", bool raw=false);
//----------------------------------------------------------------------------//

    /**
     * Enables the dEdx Interpolation for every crosssection in vector
     * crosssections_
     */
    void EnableDEdxInterpolation(std::string path ="", bool raw=false);

//----------------------------------------------------------------------------//
    /**
     * Enables the dNdx Interpolation for every crosssection in vector
     * crosssections_
     */
    void EnableDNdxInterpolation(std::string path ="", bool raw=false);

//----------------------------------------------------------------------------//
    /**
     * Disables the dEdx Interpolation for every crosssection in vector
     * crosssections_
     */
    void DisableDEdxInterpolation();

//----------------------------------------------------------------------------//
    /**
     * Disables the dNdx Interpolation for every crosssection in vector
     * crosssections_
     */
    void DisableDNdxInterpolation();

//----------------------------------------------------------------------------//
    /**
     * Disables the Interpolation including dEdx and dNdx for
     * every crosssection in vector crosssections_
     */
    void DisableInterpolation();

//----------------------------------------------------------------------------//
    /**
     * Disables the Interpolation for calculation the exact particle time
     */
    void DisableParticleTimeInterpolation();
//----------------------------------------------------------------------------//

    /**
     * Enables the Interpolation for calculation the exact particle time
     */
    void EnableParticleTimeInterpolation(std::string path="", bool raw=false);
//----------------------------------------------------------------------------//

    void EnableLpmEffect();

//----------------------------------------------------------------------------//

    void DisableLpmEffect();

//----------------------------------------------------------------------------//

    void EnableContinuousRandomization();

//----------------------------------------------------------------------------//

    void DisableContinuousRandomization();

//----------------------------------------------------------------------------//

    void EnableScattering();

//----------------------------------------------------------------------------//

    void DisableScattering();

//----------------------------------------------------------------------------//

    void EnableExactTimeCalculation();

//----------------------------------------------------------------------------//

    void DisableExactTimeCalculation();

//----------------------------------------------------------------------------//
    /*!
    * function for range calculation for given energy - interface to Integral;
    * \f[f(E) =- \frac{1}{ \frac{dE}{dx}\big|_{Ioniz} +\frac{dE}{dx}\big|
    * _{Brems}+\frac{dE}{dx}\big|_{PhotoNuc}+\frac{dE}{dx}\big|_{epair}  }\f]
    * \return range calculation for given energy [cm/E]
    * \param E energy [MeV]
    */

    double FunctionToIntegral(double energy);

//----------------------------------------------------------------------------//

    // destructors

    ///@brief Crush this CrossSections.
    virtual ~ProcessCollection();

//----------------------------------------------------------------------------//
    //Getter

    Decay* GetDecay() const {
        return decay_;
    }

    int GetLocation() const {
        return location_;
    }

    std::vector<CrossSections*> GetCrosssections() const {
        return crosssections_;
    }

    EnergyCutSettings* GetCutSettings() const {
        return cut_settings_;
    }

    bool GetDebug() const {
        return debug_;
    }

    bool GetDoInterpolation() const {
        return do_interpolation_;
    }

    bool GetDoRandomization() const {
        return do_continuous_randomization_;
    }

    bool GetDoScattering() const {
        return do_scattering_;
    }

    double GetIni() const {
        return ini_;
    }

    Integral* GetIntegral() const {
        return integral_;
    }

    Interpolant* GetInterpolant() const {
        return interpolant_;
    }

    Interpolant* GetInterpolantDiff() const {
        return interpolant_diff_;
    }

    bool GetLpmEffectEnabled() const {
        return lpm_effect_enabled_;
    }

    ContinuousRandomization* GetContinuousRandomization() const {
        return randomizer_;
    }

    Scattering* GetScattering() const {
        return scattering_;
    }

    bool GetEnableRandomization() const {
        return enable_randomization_;
    }

    Medium* GetMedium() const {
        return medium_;
    }

    int GetOrderOfInterpolation() const {
        return order_of_interpolation_;
    }

    PROPOSALParticle* GetParticle() const {
        return particle_;
    }

    Geometry* GetGeometry() const {
        return geometry_;
    }

    double GetDensityCorrection() const {
        return density_correction_;
    }
//----------------------------------------------------------------------------//
    //Setter

    void SetCrosssections(std::vector<CrossSections*> crosssections);
    void SetCutSettings(EnergyCutSettings* cutSettings);
    void SetDebug(bool debug);
    void SetDoInterpolation(bool doInterpolation);
    void SetIni(double ini);
    void SetIntegral(Integral* integral);
    void SetInterpolant(Interpolant* interpolant);
    void SetInterpolantDiff(Interpolant* interpolantDiff);
    void SetLpmEffectEnabled(bool lpmEffectEnabled);
    void SetMedium(Medium* medium);
    void SetOrderOfInterpolation(int orderOfInterpolation);
    void SetParticle(PROPOSALParticle* particle);
    void SetLocation(int location);
    void SetGeometry(Geometry* geometry);
    void SetDensityCorrection(double density_correction);
    void SetEnableRandomization(bool enable_randomization);

    PROPOSALParticle *GetBackup_particle() const;
    void SetBackup_particle(PROPOSALParticle *backup_particle);
    void RestoreBackup_particle();
};

}


#endif //ProcessCollection_H
