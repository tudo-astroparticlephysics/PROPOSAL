#ifndef EXPCORR_H
#define EXPCORR_H

#include "generation/gen/CosCorr.h"
#include "PROPOSAL/Integral.h"
#include <string>
#include "PROPOSAL/PhysicsModel.h"



/**
 * This class calculates atmosphere-dependent quantities for a given average first interaction depth.
 */

class ExpCorr :public PhysicsModel{



private:
    int func;
    double xs;
    double Norm;
    Integral *I;

public:
    double h0;
    double X0;

    CosCorr *C;


    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize class with atmospheric model and ground elevation z0 in [km].
     * h0 is the average production height in km or in g/cm^2 if negative.
     */

    ExpCorr(int model, double h0, double z0);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Prints zenith angle profile of production height, cos*, total X, track length.
     */

   // static void main(std::string[] args);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Common function for calculation of atmospheric averages:
     * <ol>
     * <li> func=1: muon production height
     * <li> func=2: 1/cos~1/rho
     * <li> func=3: muon track length
     * <li> all others: normalization.
     * </ol>
     */

    double getExp(int func, double x);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Common function for calculation of atmospheric averages - interface to Integral.
     */

    double function(double X);

};


#endif // EXPCORR_H
