#ifndef INTFLUX_H
#define INTFLUX_H

#include "PROPOSAL/PhysicsModel.h"
#include "generation/gen/EarthModel.h"
#include <string>
#include "generation/gen/ExpCorr.h"
#include "PROPOSAL/Integral.h"
#include "generation/gen/NeutrinoTot.h"
#include "PROPOSAL/Interpolate.h"
#include <vector>

/**
 * Class defines Gaisser-like muon and muon- and electron-neutrino energy spectra with cos*, muon energy loss, and decay corrections.
 */

class IntFlux :public PhysicsModel{

private:
    double z0;
    double X0;
    double R;
    int P;
    int M;
    int S;

    double xpar[4][5];
//    ={
//        /* Standard US Atmosphere average parameters */
//        { 0.00851164, 0.124534, 0.059761, 2.32876, 19.279 },
//        { 0.102573, -0.068287, 0.958633, 0.0407253, 0.817285 },
//        { -0.017326, 0.114236, 1.15043, 0.0200854, 1.16714 },
//        { 1.3144, 50.2813, 1.33545, 0.252313, 41.0344 }
//    };


    double fpar[9][5];
//    ={
//        {0.701, 2.715, 1, 1, 1},
//        {0.340, 2.720, 1, 1, 1},
//        {0.367, 2.712, 1, 1, 1},
//        {0.646, 2.684, 1, 1, 1},
//        {0.352, 2.680, 1, 1, 1},
//        {0.310, 2.696, 1, 1, 1},
//        {0.828, 2.710, 1, 1, 1},
//        {0.465, 2.719, 1, 1, 1},
//        {0.472, 2.741, 1, 1, 1}
//        // { 0.701459, 2.71461, 1, 1, 1 }; all corrections
//        // { 0.590968, 2.69506, 1, 0, 0 }; just cos* correction
//        // { 0.594662, 2.69671, 0, 0, 0 }; original formula
//    };

    double dG;
    double xx, xs, xo, xl, xn;
    std::string name;
    ExpCorr *eC;
    Integral *I;

    int hcol;

    void init(int type, int model, int flag, double h0, double z0);

    //----------------------------------------------------------------------------------------------------//

    void flSet(double x);

    //----------------------------------------------------------------------------------------------------//

    double h(double x);

    //----------------------------------------------------------------------------------------------------//

    double c(double x);

    //----------------------------------------------------------------------------------------------------//

    double o(double x);

    //----------------------------------------------------------------------------------------------------//

    double l(double x);

    //----------------------------------------------------------------------------------------------------//



public:
    double D;
    double Rc;

    static int sM;
    static double Efr;
    static double mF;
    static double gD;
    static double gEcut;
    static EarthModel *eM;
    static NeutrinoTot *nT;


    Interpolate *J;
    bool jt;


    static std::vector<Interpolate*> Je;
    static bool je;



    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize class with particle type (all mu/mu-/mu+/all nu_mu/nu_mu/~nu_mu/all nu_e/nu_e/~nu_e/),
     * model and ground elevation z0 in [km]. flag switches between two cos* calculation algorithms.
     * h0 is the average production height in km or in g/cm^2 if negative.
     */

    IntFlux(int type, int model, int flag, double h0, double z0);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize the class for spectrum index correction dG and prompt to pion muon component ratio Rc.
     */

    IntFlux(int type, int model, int flag, double h0, double z0, double dG, double Rc, double D);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Parametrize h(x) c(x) o(x) l(x).
     */

    void interpolate(std::string name);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates 3 example distrbutions.
     */

   // static void main(std::string[] args);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates integral flux above energy E in [GeV] at x=cos(zenith angle).
     */

    double getIntFlux(double E, double x);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates particle energy above threshold E0 in [GeV] at x=cos(zenith angle) as a function of a random number rnd.
     */

    double getE(double E0, double x, double rnd);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates differential flux for energy E in [GeV] at x=cos(zenith angle).
     */

    double getfl(double E, double x);



    /**
     * Energy distribution as a function of energy x in [GeV]; zenith angle is set in a private variable - interface to Integral.
     */

    double function(double x);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Parametrizes the integral of this class.
     */

    void interpolate(double E0);

    //----------------------------------------------------------------------------------------------------//

    /**
     * 2d parametrization - interface to Interpolate
     */

    double functionInt(double x, double E);

    //----------------------------------------------------------------------------------------------------//



    /**
     * 1d parametrization - interface to Interpolate
     */

    double functionInt(double x);

};


#endif // INTFLUX_H
