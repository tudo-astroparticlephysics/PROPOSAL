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



/*! \class ProcessCollection ProcessCollection.h "CrossSections.h"
    \brief initializes all cross sections and keeps references to them
 */

class ProcessCollection
{

protected:

    int         order_of_interpolation_;
    bool        do_interpolation_;
    bool        lpm_effect_enabled_;
    double      ini_;
    bool        debug_;

    Interpolant*        interpolant_;
    Interpolant*        interpolant_diff_;
    Particle*           particle_;
    Medium*             medium_;
    Integral*           integral_;
    EnergyCutSettings*  cut_settings_;

    std::vector<CrossSections*> crosssections_;



public:



//----------------------------------------------------------------------------//

    // constructors

    ProcessCollection();
    ProcessCollection(const ProcessCollection&);
    ProcessCollection& operator=(const ProcessCollection&);
//----------------------------------------------------------------------------//

    /// @brief  initializes all cross sections,

    ProcessCollection(Particle *particle, Medium *medium, EnergyCutSettings* cut_settings);

//----------------------------------------------------------------------------//

    // Memberfunctions

    /*!
    * function for range calculation for given energy - interface to Integral;
    * \f[f(E) =- \frac{1}{ \frac{dE}{dx}\big|_{Ioniz} +\frac{dE}{dx}\big|
    * _{Brems}+\frac{dE}{dx}\big|_{PhotoNuc}+\frac{dE}{dx}\big|_{epair}  }\f]
    * \return range calculation for given energy [cm/E]
    * \param E energy [MeV]
    */

    double FunctionToIntegral(double energy);

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
    double GetDx(double ei, double ef, double dist);

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

    double GetEf(double ei, double dist);

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
    void EnableInterpolation();

//----------------------------------------------------------------------------//

    void DisableInterpolation();

//----------------------------------------------------------------------------//

    void EnableLpmEffect();

//----------------------------------------------------------------------------//

    void DisableLpmEffect();

//----------------------------------------------------------------------------//

    // destructors

    ///@brief Crush this CrossSections.
    virtual ~ProcessCollection();

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
