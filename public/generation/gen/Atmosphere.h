#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H

#include <string>
#include "PROPOSAL/PhysicsModel.h"

/**
 * Atmosphere parameter implementation.
 */

class Atmosphere :public PhysicsModel{

private:
    int w;
    int S;
    double z0;

    double h1;
    double h2;

    double A[3][4];
//    ={
//        {-186.5562, -94.919, 0.61289, 0}, /* Std US */
//        {-142.801, -70.1538, 1.14855, 0}, /* Spl Oct 01 */
//        {-163.331, -65.3713, 0.402903, 0} /* Spl Jul 01 */
//    };

    double B[3][4];
//    ={
//        {1222.6562, 1144.9069, 1305.5948, 540.1778},
//        {1177.19, 1125.11, 1304.77, 433.823},
//        {1183.70, 1108.06, 1424.02, 207.595}
//    };
    double C[3][4];
//    ={
//        {9.9418638, 8.7815355, 6.3614304, 7.7217016},
//        {8.61745, 7.65925, 5.81351, 7.75155},
//        {8.75221, 7.53213, 5.45846, 7.93043}
//    };
    double H[3][4];
//    ={
//        {4, 10, 40, 100},
//        {4, 10, 40, 100},
//        {4, 10, 40, 100}
//    };
    double X_[4];
    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize class with model and ground elevation z0 in [km].
     */

public:
    Atmosphere(int model, double z0);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Prints atmospheric profiles.
     */

   // static void main(std::string[] args);

    //----------------------------------------------------------------------------------------------------//



    //----------------------------------------------------------------------------------------------------//

    /**
     * Returns altitude in [km] at a function of mass overburden in [g/cm^2].
     */

    double h(double x);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Returns mass overburden in [g/cm^2] at a function of altitude in [km].
     */

    double X(double h);
    //----------------------------------------------------------------------------------------------------//

    /**
     * Returns density in [g/cm^2/km] at a function of mass overburden in [g/cm^2].
     */

    double dXdhX(double x);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Returns density in [g/cm^2/km] at a function of altitude in [km].
     */

    double dXdh(double h);

};

#endif // ATMOSPHERE_H
