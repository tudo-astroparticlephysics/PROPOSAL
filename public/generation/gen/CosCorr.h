#ifndef COSCORR_H
#define COSCORR_H

#include "generation/gen/EarthModel.h"
#include "generation/gen/Atmosphere.h"
#include <string>
#include "PROPOSAL/PhysicsModel.h"


/**
 * This class calculates atmosphere-dependent quantities for a fixed first interaction depth.
 */

class CosCorr :public PhysicsModel{

private:
    double R;
    Integral *I;


    double st;
    double sum;

public:

    Atmosphere *A;
    double h0;
    double X0;
    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize class with atmospheric model and ground elevation z0 in [km].
     * h0 is the average production height in km or in g/cm^2 if negative.
     */

    CosCorr(int model, double h0, double z0);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Prints zenith angle profile of production height, cos*, total X, track length.
     */

   // static void main(std::string[] args);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates altitude in [km] corresponding to depth X0 in [g/cm^2] as a function of x=cos(zenith angle).
     */

    double geth(double x, double X0);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates atmospheric depth in [g/cm^2] as a function of x=cos(zenith angle).
     */

    double getX(double x);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates atmospheric depth in [g/cm^2] above altitude h in [km] as a function of x=cos(zenith angle).
     */

    double getX(double x, double h);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates atmospheric density in [g/cm^2/km] above altitude h in [km] as a function of x=cos(zenith angle).
     */

    double dXdh(double x, double h);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates path length in [km] to altitude h in [km] as a function of x=cos(zenith angle).
     */

    double getx(double x, double h);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates derivative of mass overburden wrt path distance in [km] as a function of x=cos(zenith angle).
     */

    double function(double x);
};

#endif // COSCORR_H
