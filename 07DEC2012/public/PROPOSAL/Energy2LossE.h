/*
 * Energy2LossE.h
 *
 *  Created on: 24.06.2010
 *      Author: koehne
 */

#ifndef ENERGY2LOSSE_H_
#define ENERGY2LOSSE_H_

#include "PROPOSAL/Energy2Loss.h"
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Interpolate.h"


/**
 * class contains functions for evaluation of the
 * spread of the continuous energy losses
 */

class Energy2LossE: public  Energy2Loss
{

protected:
    Integral *integral_;
    bool df;
    bool jt;


//----------------------------------------------------------------------------//



public:

    Interpolate *interpolateJ_;
    Interpolate *interpolateJdf_;

    /**
    * initializes an integral class separate from that in Propagate
    */
//----------------------------------------------------------------------------//
    //Constructor

    Energy2LossE(Energy2Loss *e2loss);

//----------------------------------------------------------------------------//
    //Memberfunctions

    /**
    * total energy^2 losses - interface to Integral
    */

    double function(double E);

//----------------------------------------------------------------------------//

    /**
    * contribution of everything to -dE2/de
    */

    double dE2de(double ei, double ef);

//----------------------------------------------------------------------------//

    /**
    * 1d parametrization - interface to Interpolate
    */

    double functionInt(double e);

//----------------------------------------------------------------------------//
    // Getter

    Interpolate* get_interpolateJ()
    {
        return interpolateJ_;
    }

    Interpolate* get_interpolateJdf()
    {
        return interpolateJdf_;
    }

//----------------------------------------------------------------------------//
    // Setter

    void set_jt(bool newJT);
    void set_df(bool newDF);

};



#endif /* ENERGY2LOSSE_H_ */
