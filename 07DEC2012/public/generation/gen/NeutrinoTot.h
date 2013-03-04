#ifndef NEUTRINOTOT_H
#define NEUTRINOTOT_H

#include "PROPOSAL/PhysicsModel.h"
#include "PROPOSAL/Interpolate.h"
#include "PROPOSAL/Integral.h"
#include "generation/gen/NeutrinoInt.h"
#include <string>
#include <vector>

/**
 * Class contains functions for calculation of neutrino interaction cross sections.
 */

class NeutrinoTot : public PhysicsModel{

private:
    NeutrinoInt *N;
    std::vector<Integral*> I;
    double *H;

    double E;
    bool cc, nu;
    double rnd, rns;

public:
    std::vector<Interpolate*> J;
    bool jt;
    std::vector<Interpolate*> Jo;
    //----------------------------------------------------------------------------------------------------//

    /**
     * Class constructor.
     */

    NeutrinoTot();

    //----------------------------------------------------------------------------------------------------//

    /**
     * Charged and Neutral current interaction for neutrinos and antineutrinos
     */

    double dSdy(double E, bool cc, bool nu);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Charged and Neutral current interaction for neutrinos and antineutrinos
     */

    double dSdy(double E, bool cc, bool nu, double rnd);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Returns the energy transferred to the hardonic state
     */

    double e(bool cc, bool nu);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Glashow resonance cross section
     */

    double dSdy(double E);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Glashow resonance cross section
     */

    double dSdy(double E, double rnd);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Glashow resonance energy lost by neutrino
     */

    double e(double m);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Neutrino oscillation probability (mu->tau); units: E [GeV] L [m]
     */

    double Pmt(double E, double L);

    //----------------------------------------------------------------------------------------------------//

    /**
     * returns the neutrino cross section dSdy - interface to Integral.
     */

    double function(double y);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Prints the cross sections in [cm^2] as a function of energy in [GeV].
     */

   // static void main(std::string[] args);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Parametrizes the integral of this class.
     */

    void interpolate(std::string name);

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

    //----------------------------------------------------------------------------------------------------//



    /**
     * 1d parametrization - interface to Interpolate
     */

    double functionInt(double E);
};


#endif // NEUTRINOTOT_H
