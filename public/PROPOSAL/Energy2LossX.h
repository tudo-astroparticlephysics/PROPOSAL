/*
 * Energy2LossX.h
 *
 *  Created on: 24.06.2010
 *      Author: koehne
 */

#ifndef ENERGY2LOSSX_H_
#define ENERGY2LOSSX_H_

#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Energy2Loss.h"

/**
 * class contains functions for evaluation
 * of the spread of the continuous energy losses
 */

class Energy2LossX : public Energy2Loss
{

protected:
    Integral *integral_;
    int name;
    bool jt_;


//----------------------------------------------------------------------------//

public:
    Interpolate *interpolateJ_;

//----------------------------------------------------------------------------//
    /**
    * initializes an integral class separate from that in Propagate
    */

    Energy2LossX(Energy2Loss e2loss);

//----------------------------------------------------------------------------//

    /**
    * total energy^2 losses - interface to Integral
    */


    double function(double v);

//----------------------------------------------------------------------------//

    /**
    * contribution of everything to -dE2/dx
    */

    double dE2dx();

//----------------------------------------------------------------------------//


    /**
    * 1d parametrization - interface to Interpolate
    */

    double functionInt(double e);

//----------------------------------------------------------------------------//
    // Getter

    bool get_jt()
    {
        return jt_;
    }

    Interpolate get_interpolateJ()
    {
        return *interpolateJ_;
    }

//----------------------------------------------------------------------------//
    // Setter

    void set_jt(bool newJT);
};


#endif /* ENERGY2LOSSX_H_ */
