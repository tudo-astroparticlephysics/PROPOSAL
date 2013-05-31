/*
 * ProcessCollection.h
 *
 *  Created on: 29.04.2013
 *      Author: koehne
 */

#ifndef ProcessCollection_H
#define ProcessCollection_H

#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Particle.h"
#include "PROPOSAL/CrossSections.h"
#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Epairproduction.h"
#include "PROPOSAL/Ionization.h"
#include "PROPOSAL/Photonuclear.h"
#include "PROPOSAL/Decay.h"
#include "PROPOSAL/Medium.h"
#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/Interpolant.h"
#include <vector>
#include <utility>
#include "PROPOSAL/ContinuousRandomization.h"



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

    bool        do_continuous_randomization; //!< Enables the randomization of continuous energy losses

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
    Particle*           particle_;
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
//----------------------------------------------------------------------------//

    /// @brief  initializes all cross sections,

    ProcessCollection(Particle *particle, Medium *medium, EnergyCutSettings* cut_settings);

//----------------------------------------------------------------------------//

    // Memberfunctions
    void swap(ProcessCollection &collection);

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
    std::pair<double,std::string> MakeStochasticLoss();

//----------------------------------------------------------------------------//

    /**
     *  Makes Decay
     *
     *  \return pair of product energy [MeV] and kind of product
     */
    std::pair<double,std::string> MakeDecay();

//----------------------------------------------------------------------------//
    /**
     * Enables the Interpolation including dEdx and dNdx for
     * every crosssection in vector crosssections_
     */
    void EnableInterpolation();
//----------------------------------------------------------------------------//

    /**
     * Enables the dEdx Interpolation for every crosssection in vector
     * crosssections_
     */
    void EnableDEdxInterpolation();

//----------------------------------------------------------------------------//
    /**
     * Enables the dNdx Interpolation for every crosssection in vector
     * crosssections_
     */
    void EnableDNdxInterpolation();

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

    void EnableLpmEffect();

//----------------------------------------------------------------------------//

    void DisableLpmEffect();

//----------------------------------------------------------------------------//

    void EnableContinuousRandomization();

//----------------------------------------------------------------------------//

    void DisableContinuousRandomization();

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
        return do_continuous_randomization;
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

	Medium* GetMedium() const {
		return medium_;
	}

	int GetOrderOfInterpolation() const {
		return order_of_interpolation_;
	}

	Particle* GetParticle() const {
		return particle_;
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
	void SetParticle(Particle* particle);
};



#endif //ProcessCollection_H
