#ifndef NEUTRINOINT_H
#define NEUTRINOINT_H

#include "PROPOSAL/Medium.h"
#include "generation/gen/CteqPDF.h"
#include <string>
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Interpolate.h"
#include <vector>

/**
 * Class contains functions for calculation of neutrino interaction differential cross sections.
 */

class NeutrinoInt : public PhysicsModel{



private:
    CteqPDF *F;
    double M, d1, d2, d3;

    std::string name;
    double y, E;
    bool cc, nu;
    Integral *I;


    //----------------------------------------------------------------------------------------------------//

    /**
     * returns the neutrino (nu=true) and anti neutrino (nu=false) charged (cc=true) and neutral (cc=false) current
     * interaction cross sections as functions of q(x,y) and y and energy E in [GeV].
     */



public:
    double dS2dxdy(double q, double y, double E, bool cc, bool nu);

    Medium *m;

    double RA;
    double RZ;


    std::vector<Interpolate*> J;
    bool jt;

    //----------------------------------------------------------------------------------------------------//

    /**
     * Class constructor. Sets up the values of cross section parameters and calls the CteqPDF constructor.
     */

    NeutrinoInt();


    //----------------------------------------------------------------------------------------------------//

    /**
     * returns the neutrino cross sections integrated over x as functions of y and energy E in [GeV].
     */

    double dSdy(double y, double E, bool cc, bool nu);

    //----------------------------------------------------------------------------------------------------//

    /**
     * returns the neutrino cross section dS2dxdy - interface to Integral.
     */

    double function(double x);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Parametrizes the integral of this class.
     */

    void interpolate();

    //----------------------------------------------------------------------------------------------------//

    /**
     * 2d parametrization - interface to Interpolate
     */

    double functionInt(double E, double y);

};

#endif // NEUTRINOINT_H
