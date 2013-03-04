#ifndef ATMFLUX_H
#define ATMFLUX_H

#include "PROPOSAL/PhysicsModel.h"
#include "generation/gen/EarthModel.h"
#include "generation/gen/NeutrinoTot.h"
#include "generation/gen/IntFlux.h"
#include "PROPOSAL/Amanda.h"
#include <vector>
#include <string>
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/PROPOSALParticle.h"
#include <deque>

/**
 * This is the main lepton generator class.
 */

class AtmFlux :public PhysicsModel{


private:
    double R;
    double *tf, *af, *sF, tF;
    double nF, mF, xn;
    double tT;
    double t0;
    long evt;
    bool OSC;
    bool GEN;
    bool PROP;

    double R0;

    int gens, lept;
    double px, py, pz, th, phi;


    bool decayflag;
    std::vector<IntFlux*> iF;
    NeutrinoTot *nT;
    Amanda *mP, *tP;
    Integral *I;
   // Random rand;

    std::vector<PROPOSALParticle*> I3hist;
    bool f2k;
    double *Ecut, *A, *G, *a, *g;
    std::string tdir;
    std::string dtct;
    std::string name[13];

    //----------------------------------------------------------------------------------------------------//

    /**
     * general-purpose lepton propagator.
     */

    void propagate(std::string type, int igen, int gens, double x, double y, double z, double th, double phi, double r, double E, double t);

    //----------------------------------------------------------------------------------------------------//

    /**
     * general-purpose lepton propagator.
     */

    void prop(int lept, int igen, double E, double x, double y, double z, double t, double nf);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Constructs the lepton from the given information.
     */

    void putOut(int lept, int igen, double x, double y, double z, double t, double l, double E, PROPOSALParticle* p);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Copies the lepton into the output stream.
     */

    void lptOut(PROPOSALParticle* p);

    //----------------------------------------------------------------------------------------------------//

    double getTotFlux(double x);

    //----------------------------------------------------------------------------------------------------//

    double getTotFlux(double x, double rnd);


public:
    double elost, eini, thini;
    std::string Name;
    EarthModel *eM;
    double D;
    double length, radius;


    std::vector<Interpolate*> J;
    bool jt;

    AtmFlux();
    //----------------------------------------------------------------------------------------------------//

    /**
     * Class initializer, command-line option parser. Call with "-help" to list all options.
     */

    long initAtmFlux(std::deque<std::string> args, bool f2k);

    //----------------------------------------------------------------------------------------------------//

    /**
     * begins the f2k stream.
     */

    void beginAtmFlux();

    //----------------------------------------------------------------------------------------------------//

    /**
     * concludes the f2k stream.
     */

    void endAtmFlux();

    //----------------------------------------------------------------------------------------------------//

    /**
     * F2k stream lepton generator.
     */

   // static void main(std::string[] args);

    //----------------------------------------------------------------------------------------------------//

    /**
     * I3m initializer.
     */

    void setup(std::deque<std::string> args);

    //----------------------------------------------------------------------------------------------------//

    /**
     * I3m initializer.
     */

    void setup(std::string args);
    //----------------------------------------------------------------------------------------------------//

    /**
     * I3m function - creates particle set of the next event.
     */

    std::vector<PROPOSALParticle*> createNext();

    //----------------------------------------------------------------------------------------------------//

    /**
     * creates the next event.
     */

    void findNext();
    //----------------------------------------------------------------------------------------------------//

    /**
     * general-purpose lepton propagator.
     */

    std::vector<PROPOSALParticle*> propagate(PROPOSALParticle* p);
    //----------------------------------------------------------------------------------------------------//

    /**
     * Main parser for the f2k file streams.
     */

    void mmcf2k();




    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates if the coordinates are inside the detector
     */

    bool insDet(double x, double y, double z);



    //----------------------------------------------------------------------------------------------------//

    /**
     * Function for flux integration over zenith angles - interface to Integral.
     */

    double function(double x);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Parametrizes the integral of this class.
     */

    void interpolate();

    //----------------------------------------------------------------------------------------------------//


    /**
     * 1d parametrization - interface to Interpolate
     */

    double functionInt(double x);

};

#endif // ATMFLUX_H
