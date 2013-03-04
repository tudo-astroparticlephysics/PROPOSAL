/*! \file   PROPOSALParticle.h
*   \brief  Header file for the PROPOSALParticle routines.
*
*   For more details see the class documentation.
*
*   \date   21.06.2010
*   \author Jan-Hendrik Köhne
*/


#ifndef PROPOSALParticle_H
#define PROPOSALParticle_H
#include "PROPOSAL/Propagate.h"
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Interpolate.h"
#include "PROPOSAL/Output.h"
#include <string>
#include "PROPOSAL/PhysicsModel.h"
#include "PROPOSAL/Scattering.h"
#include "vector"


class Integral;
class Interpolate;
class Scattering;

/**
  * \brief This class provides the main particle properties and functions.
  *
  * All coordinates, angles and physical values are stored in this class.
  * Functions to calculate the time- and space differences are implemented.
  */


class PROPOSALParticle : public PhysicsModel
{


public:

    double r;                   //!< propagation distance [cm]   - Set to 0 in Constructor
    double x;                   //!< x-coordinate [cm]           - Set to 0 in Constructor
    double y;                   //!< y-coordinate [cm]           - Set to 0 in Constructor
    double z;                   //!< z-coordinate [cm]           - Set to 0 in Constructor
    double t;                   //!< age [sec]                   - Set to 0 in Constructor
    double theta;               //!< zenith of the momentum in [deg]    - Set to 0 in Constructor
    double phi;                 //!< azimuth of the momentum in [deg]   - Set to 0 in Constructor

    long double costh;          //!< cos(theta)     - Set to 1 in Constructor
    long double sinth;          //!< sin(theta)     - Set to 0 in Constructor
    long double cosph;          //!< cos(phi)       - Set to 1 in Constructor
    long double sinph;          //!< sin(phi)       - Set to 1 in Constructor

    double p;                   //!< momentum [MeV]             - Set to 0 in Constructor
    double p2;                  //!< momentum square [MeV]      - Set to 0 in Constructor
    double e;                   //!< energy [MeV]               - Set to 0 in Constructor
    double m;                 	//!< mass [MeV]         - Set to Muon Mass in Constructor:      m = MMU
    double l;                   //!< lifetime [sec]     - Set to Muon Liftime in Constructor:   l = LMU
    double c;                   //!< charge             - Set to Muon charge in Constructor:    c = 1

    std::string name;          	//!< name of the particle - Presetted to "mu"
    double low;                	//!< energy below which the particle is lost [MeV]    - Set to m in Constructor
    int type;                  	//!< particle type: 1 for muon, 2 for tau, 3 for electron - default is 1 (made in Constructor)

    int igen;                   //!< parent particle id         - Set to 0 in Constructor
    int gens;                   //!< particle id                - Set to 1 in Constructor

    Scattering  *scattering_;   //!< scattering class
    Propagate   *propagate_;    //!< propagate class
    Integral    *integral_;     //!< integral class

    double xi;                  //!< x-coordinate entry Point [m]           - Set to 0 in Constructor
    double yi;                  //!< y-coordinate entry Point [m]           - Set to 0 in Constructor
    double zi;                  //!< z-coordinate entry Point [m]           - Set to 0 in Constructor
    double ti;                  //!< t-coordinate entry Point [sec]         - Set to 0 in Constructor
    double Ei;                  //!< energy at entry point [GeV]            - Set to 0 in Constructor

    double xf;                  //!< x-coordinate exit Point [m]           - Set to 0 in Constructor
    double yf;                  //!< y-coordinate exit Point [m]           - Set to 0 in Constructor
    double zf;                  //!< z-coordinate exit Point [m]           - Set to 0 in Constructor
    double tf;                  //!< t-coordinate exit Point [sec]         - Set to 0 in Constructor
    double Ef;                  //!< energy at exit point [GeV]            - Set to 0 in Constructor

    double xc;                  //!< x-coordinate at point of closest approach [m]           - Set to 0 in Constructor
    double yc;                  //!< y-coordinate at point of closest approach [m]           - Set to 0 in Constructor
    double zc;                  //!< z-coordinate at point of closest approach [m]           - Set to 0 in Constructor
    double tc;                  //!< t-coordinate at point of closest approach [sec]         - Set to 0 in Constructor
    double Ec;                  //!< energy at at point of closest approach [GeV]            - Set to 0 in Constructor

    double Elost;              	//!< energy lost in the detector volume [GeV] - Preset to 0 in Constructor

    bool df;
    Interpolate *interpolateJ_;     //!< interpolate object for integrated function
    Interpolate *interpolateJdf_;   //!< interpolate object for differential function

    bool    jt;    //!< interpolate values?

    Output  *output;      //!< output object


//----------------------------------------------------------------------------//

    /**
     * \brief Default Constructor
     *
     * Constructor which sets "default" settings.
     */
    PROPOSALParticle();
//----------------------------------------------------------------------------//

    /*!
     * initialize particle
     *
     * \param   *pr     Propagate object which will be set
     * \param   name    particle type e.g. "mu"
     */
    PROPOSALParticle(Propagate *pr, std::string name);

//----------------------------------------------------------------------------//
    /*!
     * \brief Create particle from another particle.
     *
     * This function ist mostly used to store particle information.
     *
     * \param igen      parent particle id
     * \param gens      particle id
     * \param name      particle type
     * \param x         x-coordinate
     * \param y         y-coordinate
     * \param z         z-coordinate
     * \param theta     theta angle
     * \param phi       phi angle
     * \param e         particle energy
     * \param t         particle time
     * \param r         flight distance
     * \param *p        source particle
     */

    PROPOSALParticle(int igen, int gens, std::string name, double x, double y, double z, double theta, double phi, double e, double t, double r, PROPOSALParticle *p);

//----------------------------------------------------------------------------//
    /*!
     * \brief Create particle with properties
     *
     * This function ist mostly used to store particle information.
     *
     * \param igen      parent particle id
     * \param gens      particle id
     * \param name      particle type
     * \param x         x-coordinate
     * \param y         y-coordinate
     * \param z         z-coordinate
     * \param theta     theta angle
     * \param phi       phi angle
     * \param e         particle energy
     * \param t         particle time
     * \param r         flight distance
     */
    PROPOSALParticle(int igen, int gens, std::string name, double x, double y, double z, double theta, double phi, double e, double t, double r);

//----------------------------------------------------------------------------//
    /*!
     * \brief Create particle with properties
     *
     * This function ist mostly used to store particle information.
     *
     * \param aname     particle type
     * \param x         x-coordinate
     * \param y         y-coordinate
     * \param z         z-coordinate
     * \param theta     theta angle
     * \param phi       phi angle
     * \param e         particle energy
     * \param t         particle time
     */
    PROPOSALParticle(std::string aname, double x, double y, double z, double theta, double phi, double e, double t);

//----------------------------------------------------------------------------//

    /*!
     * initialize particle by its name
     *
     * \param name      particle type
     */

    void initByName(std::string name);


//----------------------------------------------------------------------------//
    // Memberfunctions

    /*!
     * initialize the location and direction of the particle,\n
     * time in sec, x, y, z in cm, theta and phi in deg
     *
     * \param name      particle type
     * \param time      particle time
     * \param x         x-coordinate
     * \param y         y-coordinate
     * \param z         z-coordinate
     * \param theta     theta angle
     * \param phi       phi angle
     */

    void location(std::string name, double time, double x, double y, double z, double theta, double phi);

//----------------------------------------------------------------------------//
    /*!
    * advances the particle by the given distance
    *
    * \param    dr  flight distance
    * \param    ei  initial energy
    * \param    ef  final energy
    */

    void advance(double dr, double ei, double ef);

//----------------------------------------------------------------------------//
    /*!
    * sets the energy of the particle
    *
    * \param    e   particle energy
    */

    void setEnergy(double e);

//----------------------------------------------------------------------------//
    /*!
    * function for time delta calculation - interface to Integral
    *
    * \return   ???
    */

    double function(double E);

//----------------------------------------------------------------------------//
    /*!
    * time delta, corresponding to the given propagation distance
    *
    * \param    ei  initial energy
    * \param    ef  final energy
    * \return   time delta
    */

    double getdt(double ei, double ef);

//----------------------------------------------------------------------------//
    /*!
    * 1d parametrization - interface to Interpolate
    */

    double functionInt(double e);
//----------------------------------------------------------------------------//
    /*!
     * \brief This calculates the alpha with running couplings.
     *
     * \f[\alpha(Q^2) = \frac{\alpha_0}{1-\frac{2\alpha_0}{3\pi}\ln(\frac{Q}{m_e})}\f]
     * with v = Q
     *
     * \param   v   Q
     * \param   alpha
     */

    double alpha(double v);

//----------------------------------------------------------------------------//
    // Getter



    double get_energy()
    {
        return e;
    }

    double get_mass()
    {
        return m;
    }

    double get_momentum()
    {
        return p;
    }

    double get_charge()
    {
        return c;
    }

    double get_low()
    {
        return low;
    }

    Scattering *get_scattering()
    {
        return scattering_;
    }

//----------------------------------------------------------------------------//
    // destructors

    ///@brief Crush this PROPOSALParticle.
    ~PROPOSALParticle() {}


};



#endif //PROPOSALParticle_H
