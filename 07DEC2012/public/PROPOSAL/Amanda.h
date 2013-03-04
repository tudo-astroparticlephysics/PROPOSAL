#ifndef AMANDA_H
#define AMANDA_H

#include <string>
#include "PROPOSAL/Propagate.h"
#include <sstream>




/**
 * Implements muon/tau propagation through multiple media to/through the detector.
 * Current possibilities for detector geometry are cylinder, cuboid, and sphere.
 * The origin is at the center of the detector.
 */
class Output;

class Amanda
{

    private:

        static bool lfix;           // static instead of const static?!
        static double LENGTH;       // static instead of const static?!
        static double RADIUS;       // static instead of const static?!
        static double WIDTH;        // static instead of const static?!
        static double HEIGHT;       // static instead of const static?!
        static double BIG_;          // static instead of const static?!

        bool SURF;
        bool FACE;
        bool USER;
        bool USFI;
        bool SDEC;
        bool RECC;
        double zset;
        std::string type;

        int mediamax;
        std::string mediadef;

        double vcut[2];
        double ecut[2];
        double vaux;
        double eaux;
        std::string med;
        std::string muta;
        std::string usna;
        double elow;
        double ebig;

        bool timef;
        bool conti[2];
        bool lpmef;
        bool scatt;
        bool frho;
        bool rfix;
        bool amasim;

        int bspar;
        int pncrs;
        int pncbb;
        int pncsh;

        double crsci;
        double crscb;
        double crscp;
        double crsce;
        double crscd;
        double drho;

        int romb;
        long seed;
        bool SEED;

        std::string tdir;
        std::vector<Propagate*> p1, p2, p3;
        double surf;
        double emax;
        std::stringstream hist;
        //Vector userbf;
        PROPOSALParticle *pI;
        int fnu;
        int *fnm;
        double *fnx;
        std::string prnt;
        int igen, gens, imax;
        std::string param;
        std::stringstream options;

public:
        double hx1, hx2;


        int flag;
        double e;

        //----------------------------------------------------------------------------------------------------//


        /**
         * Calculate points of intersection with the cylinder.
         */

        void setcyl(double x, double y, double z, double cosph, double sinph, double costh,double sinth, double length, double radius);


        //----------------------------------------------------------------------------------------------------//

        /**
         * Calculate points of intersection with the box.
         */

        void setbox(double x, double y, double z, double dx, double dy, double dz, double length, double width, double height);


        //----------------------------------------------------------------------------------------------------//

        /**
         * Calculate points of intersection with the sphere.
         */

        void setsph(double x, double y, double z, double cosph, double sinph, double costh, double sinth, double radius);


        //----------------------------------------------------------------------------------------------------//

        /**
         * Calculate points of intersection with the box.
         */

        void setpln(double x, double y, double z, double dx, double dy, double dz, double nx, double ny, double nz);



        //----------------------------------------------------------------------------------------------------//

        /**
         * checks the flag, corrects the units, initializes the particle and calls the propagator
         */

        double propagateTo(double h, double e, double n, int i, Propagate *p, int igen, int gens,
                                          double time, double x, double y, double z, double theta, double phi);




        /**
          * Standartconstructor
          */

        Amanda();

     std::string dtct;

     std::vector<PROPOSALParticle*> I3hist;

     bool raw;
     int gdet;
     double length, radius, width, height;
     double *rho, *srho;



     bool dw;
     double rw, hw, fw;
     long events, tracks, vertices, missed;



     int *medt;
     double *sphz;
     double *sphr;
     double *boxx;
     double *boxy;
     double *boxz;
     double *boxl;
     double *boxw;
     double *boxh;
     double *cylz;
     double *cylr;
     double *cyll;
     int medianum;
     bool DEBUG;

    //----------------------------------------------------------------------------------------------------//



    /**
     * user-block variable: z-coordinate of entry point into the detector cylinder (Z_IN)
     */

     double z1;

    /**
     * user-block variable: z-coordinate of exit point from the detector cylinder (Z_OUT)
     */

     double z2;

    /**
     * user-block variable: distance to the point of entry into the detector cylinder (R_IN)
     */

     double h1;

    /**
     * user-block variable: distance to the point of exit from the detector cylinder (R_OUT)
     */

     double h2;

    /**
     * user-block variable: x-coordinate of closest approach point (CPD_X) or -user=z point
     */

     double nx;

    /**
     * user-block variable: y-coordinate of closest approach point (CPD_Y) or -user=z point
     */

     double ny;

    /**
     * user-block variable: z-coordinate of closest approach point (CPD_Z) or -user=z point
     */

     double nz;

    /**
     * user-block variable: Energy at the entry point into the detector cylinder if positive (E_IN)
     */

     double e1;

    /**
     * user-block variable: Energy at the exit point from the detector cylinder if positive (E_OUT)
     */

     double e2;

    /**
     * user-block variable: Energy at the CPD point if positive (E_CPD)
     */

     double ec;

    /**
     * Energy lost inside the detector cylinder
     */

     double elost;



    //----------------------------------------------------------------------------------------------------//

    /**
     * Front-end for AMANDA, reads from and outputs to standard input/output in F2000 format.
     */

     void main(std::string args);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Initializes propagation for external applications, e.g.&nbsp;mmc-icetray (through jni).
     */

     void setup(std::string args);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Initializes propagation for external applications, e.g.&nbsp;mmc-icetray (through jni).
     */

     void setup(std::deque<std::string> args);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Propagates particles for the external applications. Returns an array of secondaries.
     * If "-user" option is used, fills user-block variables.
     */

     std::vector<PROPOSALParticle*> propagate(PROPOSALParticle *p);

    //----------------------------------------------------------------------------------------------------//

    /**
     * This is the command-line option parser. Call with "-help" to list all options.
     */

     bool setp(std::deque<std::string> args);


    //----------------------------------------------------------------------------------------------------//

    /**
     * History lines for the f2k file streams.
     */

     void history();


    //----------------------------------------------------------------------------------------------------//

    /**
     * Define user line block for the f2k file streams.
     */

     void defUSER();


    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize user line block for the f2k file streams.
     */

     void iniUSER();


    //----------------------------------------------------------------------------------------------------//

    /**
     * Record user line block for the f2k file streams.
     */

     void recUSER();


    //----------------------------------------------------------------------------------------------------//

    /**
     * Output user line block for the f2k file streams.
     */

     void outUSER();


    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize particle buffers for the f2k file streams.
     */

     void initBufs();

    //----------------------------------------------------------------------------------------------------//

    /**
     * Main parser for the f2k file streams.
     */
     // former string[] !!
     void mmcf2k(std::vector<std::string> args);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Main propagator routine.
     */

     double prop(int igen, int gens, double x, double y, double z, double th, double phi, double l, double e, double t);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Main propagator routine.
     */

     double prop(double x, double y, double z, double th, double phi, double l, double e, double t);




};


#endif // AMANDA_H
