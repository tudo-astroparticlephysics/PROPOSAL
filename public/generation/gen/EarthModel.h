#ifndef EARTHMODEL_H
#define EARTHMODEL_H

#include "PROPOSAL/PhysicsModel.h"
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Interpolate.h"
#include "generation/gen/Atmosphere.h"
#include <string>


/**
 * Earth density model implementation.
 */

class EarthModel :public PhysicsModel{



private:
    double z0;
    double b0;
    double D;
    double Rs;
    Atmosphere *A;
    Integral *I;

    double Xmid;

    double X0, t1, t2;

    double ts;

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates distance to the center r in [km] as a function of path length in [km]. sett must be called first.
     */


    double r(double t);
    //----------------------------------------------------------------------------------------------------//

    double jInt(double t);

    //----------------------------------------------------------------------------------------------------//

    double fInt(double X);


public:
    const static double R0=6371.3;

    int num;
    double *mRa, *mRo;
    double px, py, pz, ti, tf;
    Interpolate *J;
    bool jt;
    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize class with atmospheric model, ground elevation z0 in [km], bedrock elevation b0 in [km], and detector depth D.
     */

    EarthModel(int model, double z0, double b0, double D);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates density in [g/cm^3] as a function of distance to the center R in [km] according to the preliminary Earth model.
     */

    double rho(double R);


    //----------------------------------------------------------------------------------------------------//


    /**
     * Calculates and stores path length in [km] to the entry and exit points.
     */

    void sett(double px, double py, double pz);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates path length in [km] from the entry point. sett must be called first.
     */

    double r(double x, double y, double z);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates mass overburden in [g/cm^2] as a function of path length in [km] to the entry and exit points.
     */

    double X(double t1, double t2);

    //----------------------------------------------------------------------------------------------------//



    /**
     * Calculates mass overburden in [g/cm^2] as a function of path length in [km] to the entry and exit points.
     * Also calculates path length corresponding to the mass overburden X0 in [g/cm^2] from the entry point.
     */

    double X(double t1, double t2, double X0);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Returns the path length in [km] corresponding to the mass overburden X0 in [g/cm^2] set by a call to X(t1, t2, X0).
     */

    double t();

    //----------------------------------------------------------------------------------------------------//

    /**
     * Returns the density in [g/cm^3] as a function of the path length in [km] - interface to Integral.
     */

    double function(double x);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Returns the density in [g/cm^3] as a function of the path length in [km], gives density of the internal medium.
     */

    double intRho(double x);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Prints the density profile.
     */

//    static void main(std::string[] args);



    //----------------------------------------------------------------------------------------------------//



    /**
     * Parametrizes the integral of this class.
     */

    void interpolate();

    //----------------------------------------------------------------------------------------------------//



    /**
     * 2d parametrization - interface to Interpolate
     */

    double functionInt(double x, double t);

};


#endif // EARTHMODEL_H
