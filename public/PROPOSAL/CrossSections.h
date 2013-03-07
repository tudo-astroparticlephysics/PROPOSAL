/*
 * CrossSections.h
 *
 *  Created on: 21.06.2010
 *      Author: koehne
 */

#ifndef CrossSections_H
#define CrossSections_H

#include "PROPOSAL/Integral.h"
#include <cmath>
#include "PROPOSAL/PhysicsModel.h"
#include "PROPOSAL/PROPOSALParticle.h"
#include "PROPOSAL/Output.h"


/*! \class CrossSections CrossSections.h "CrossSections.h"
    \brief initializes all cross sections and keeps references to them
 */

class Medium;
class Interpolate;
class Decay;
class Ionizationloss;
class Epairproduction;
class Bremsstrahlung;
class Photonuclear;

class CrossSections: public PhysicsModel
{


protected:


    Bremsstrahlung *bremsstrahlung_;
    Photonuclear *photonuclear_;
    Ionizationloss *ionization_;
    Epairproduction *epairproduction_;
    Decay *decay_;
    // Mpairproduction *mpairproduction_;
    Output *o;

    int component_;		// Set to 0 in Constructor

    // those c' values are values to turn crossections on of or whatever you want to do.
    // Be careful if you change them to another value than 1

    double ci_;     // Set to 1 in Constructor
    double cb_;		// Set to 1 in Constructor
    double cp_;		// Set to 1 in Constructor
    double ce_;		// Set to 1 in Constructor
    double cd_;		// Set to 1 in Constructor
    double cm_;     // Set to 1 in Constructor

    // Value to give back upper or lower boundaries for an error value. For usual use, put them to 0

    double bremserror;
    double epairerror;
    double photoerror;
    double ionizerror;

    double e_hi;
    double e_low;
    int g;
    double ini_;
    bool df_;
    bool jt_;
    bool lpm_;



public:

    Interpolate *interpolateJ_;
    Interpolate *interpolateJdf_;
    PROPOSALParticle *particle_;
    Medium *medium_;
    Integral *integral_;
    CrossSections *cros;

//----------------------------------------------------------------------------//

    // constructors

    /// @brief  Necessary to keep class structure,
    /// since it is called from subclasses
    CrossSections();

//----------------------------------------------------------------------------//

    /// @brief  creates internal references to p and m,
    /// to be called from subclasses

    CrossSections(CrossSections *cros);

//----------------------------------------------------------------------------//

    /// @brief  initializes all cross sections,
    /// creates internal references to p and m

    CrossSections(PROPOSALParticle *p, Medium *m);

//----------------------------------------------------------------------------//

    // Memberfunctions

    /*!
    * function for range calculation for given energy - interface to Integral;
    * \f[f(E) =- \frac{1}{ \frac{dE}{dx}\big|_{Ioniz} +\frac{dE}{dx}\big|
    * _{Brems}+\frac{dE}{dx}\big|_{PhotoNuc}+\frac{dE}{dx}\big|_{epair}  }\f]
    * \return range calculation for given energy [cm/E]
    * \param E energy [MeV]
    */

    double function(double E);

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
    double getdx(double ei, double ef, double dist);

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

    double getef(double ei, double dist);

//----------------------------------------------------------------------------//

    /*!
    1d parametrization - interface to Interpolate;
    if: parameter df_ =true, \f[return=function(e)\f];
    else: \f[return= \int_{e}^{e_{low}} f(E) dE \f] with \f$
    e_{low} \f$ energy below which the particle is lost
    \param e energy [MeV]
    */
    double functionInt(double e);

//----------------------------------------------------------------------------//
    /*!
    * Sets constants needed for crosssecions. Called in constructor
    *
    */
    void SetConstants();

//----------------------------------------------------------------------------//

    // Getter

    Decay* get_decay();
    Ionizationloss* get_ionization();
    Bremsstrahlung* get_bremsstrahlung();
    Photonuclear* get_photonuclear();
    Epairproduction* get_epairproduction();

    int get_component()
    {
        return component_;
    }

    double get_ci() const
    {
        return ci_;
    }

    double get_cb() const
    {
        return cb_;
    }

    double get_cp() const
    {
        return cp_;
    }

    double get_ce() const
    {
        return ce_;
    }

    double get_cd() const
    {
        return cd_;
    }

    double get_bremserror() const
    {
        return bremserror;
    }

    double get_epairerror() const
    {
        return epairerror;
    }

    double get_photoerror() const
    {
        return photoerror;
    }

    double get_ionizerror() const
    {
        return ionizerror;
    }

    PROPOSALParticle* get_particle() const
    {
        return particle_;
    }

    Medium* get_medium() const
    {
        return medium_;
    }

    Integral* get_integral() const
    {
        return integral_;
    }

    double get_ini() const
    {
        return ini_;
    }

    bool get_df() const
    {
        return df_;
    }

    Interpolate* get_interpolateJ() const
    {
        return interpolateJ_;
    }

    Interpolate* get_interpolateJdf() const
    {
        return interpolateJdf_;
    }

    bool get_jt() const
    {
        return jt_;
    }

    bool get_lpm() const
    {
        return lpm_;
    }

//----------------------------------------------------------------------------//

    // Setter

    void set_decay(Decay *decay);
    void set_ionization(Ionizationloss *ionization);
    void set_bremsstrahlung(Bremsstrahlung *bremsstrahlung);
    void set_photonuclear(Photonuclear *photonuclear);
    void set_epairproduction(Epairproduction *epairproduction);
    void set_component(int component);
    void set_ci(double ci);
    void set_cb(double cb);
    void set_cp(double cp);
    void set_ce(double ce);
    void set_cd(double cd);
    void set_bremserror(double newbrems);
    void set_epairerror(double newepair);
    void set_photoerror(double newphoto);
    void set_ionizerror(double newioniz);
    void set_particle(PROPOSALParticle *particle);
    void set_medium(Medium *medium);
    void set_integral(Integral *integral);
    void set_ini(double ini);
    void set_df(bool df);
    void set_interpolateJ(Interpolate *interpolateJ);
    void set_interpolateJdf(Interpolate *interpolateJdf);
    void set_jt(bool jt);
    void set_lpm(bool lpm);
    
//----------------------------------------------------------------------------//

    // destructors

    ///@brief Crush this CrossSections.
    virtual ~CrossSections();


};



#endif //CrossSections_H
