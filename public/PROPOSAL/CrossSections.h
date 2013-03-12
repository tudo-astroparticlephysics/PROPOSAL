/*
 * CrossSections.h
 *
 *  Created on: 2013.03.12
 *      Author: koehne
 */

#ifndef CrossSections_H
#define CrossSections_H

/*! \class CrossSections CrossSections.h "CrossSections.h"
    \brief This is a pure virtual class
 */


class CrossSections
{


protected:


public:


//----------------------------------------------------------------------------//

    // constructors

    /// @brief  Necessary to keep class structure,
    /// since it is called from subclasses
    CrossSections();



//----------------------------------------------------------------------------//

    /// @brief  creates internal references to p and m,
    /// to be called from subclasses

    //CrossSections(CrossSections *cros);

//----------------------------------------------------------------------------//

    /// @brief  initializes all cross sections,
    /// creates internal references to p and m

    //CrossSections(PROPOSALParticle *p, Medium *m);

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



//----------------------------------------------------------------------------//

    // Setter


    
//----------------------------------------------------------------------------//

    // destructors

    ///@brief Crush this CrossSections.
    virtual ~CrossSections();


};



#endif //CrossSections_H
