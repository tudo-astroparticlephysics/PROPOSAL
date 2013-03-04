/*
 * CrossSections.h
 *
 *  Created on: 29.08.2010
 *      Author: schmitz

 */

#ifndef EPAIRSTOCHASTIC_H
#define EPAIRSTOCHASTIC_H


#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Interpolate.h"
#include <algorithm>
#include "PROPOSAL/Epairproduction.h"


/*! \class EpairStochastic EpairStochastic.h "EpairStochastic.h"
    \brief class contains functions necessary
    for calculation of e+e- pair production losses
 */

class EpairStochastic : public Epairproduction
{

protected:

    std::vector<Integral*> integral_;
    double* H_;
    double sum;
    double rnd;
    bool jt_;

//----------------------------------------------------------------------------//

public:


    Interpolate* interpolateJo_;
    Interpolate* interpolateJ_;

//----------------------------------------------------------------------------//
    //Constructor

    EpairStochastic(Epairproduction *cros);

//----------------------------------------------------------------------------//
    //Memberfunctions

    /*!
    stochastic energy losses - interface to integral;
    \f[return= c_e*e_{Pair}(v, components)\f]
    */

    double function(double v);
//----------------------------------------------------------------------------//

    /*!
    this is what the sum over all components of dN/dx is equal to;
    \f[\frac{dN}{dx} =return = c_e \cdot
    \sum_{numC} \int_{v_{Up}}^{v_{max}} e_{Pair} dv\f]
    */

    double dNdx();
//----------------------------------------------------------------------------//

    /*!
    this is what the sum over all components of dNdx is equal to;
    like dNdx(), but a double rnd is assigned, which is important for getUpper
    */

    double dNdx(double rnd);
//----------------------------------------------------------------------------//

    /*!
    this is the value Vof e=v*E, corresponding to rnd in the call to dNdx
    for the component chosen with rnd in the argument
    */

    double e(double rnd);
//----------------------------------------------------------------------------//

    /*!
    2d parametrization - interface to Interpolate;
    if \f$v_{Up}=v_{Max}\f$: \f[return=0;\f]
    else: \f$v=v_{Up}\exp\big(v\log\big(\frac{v_{Max}}{v_{Up}}\big)\big)\f$
    \f[return=\int_{v_{up}}^{v} c_e e_{Pair} dv\f]
    */


    double functionInt(double e, double v);

//----------------------------------------------------------------------------//

    double functionInt(double e);

//----------------------------------------------------------------------------//

    // Getter

    bool get_jt()
    {
        return jt_;
    }

    Interpolate* get_Jo()
    {
        return interpolateJo_;
    }

    Interpolate* get_J()
    {
        return interpolateJ_;
    }

//----------------------------------------------------------------------------//

    // Setter
    void set_jt(bool newJT)
    {
        jt_ =   newJT;
    }


};

#endif // EpairStochastic_H
