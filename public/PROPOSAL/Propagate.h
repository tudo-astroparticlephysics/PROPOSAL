/*! \file   Propagate.h
*   \brief  Header file for definition of the propagate class object.
*
*   For more details see the class documentation.
*
*   \date   05.07.2010
*   \author Jan-Hendrik Koehne
*/

#ifndef PROPAGATE_H_
#define PROPAGATE_H_

#include "PROPOSAL/PhysicsModel.h"
#include <string>
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/StandardNormal.h"
#include "PROPOSAL/PROPOSALParticle.h"

class PROPOSALParticle;
class Medium;
class CrossSections;
class Energy2Loss;
class Output;

class Interpolate;


/*! \mainpage PROPOSAL:
 * <b>PR</b>opagator with <b>O</b>ptimal <b>P</b>recision and <b>O</b>ptimized <b>S</b>peed for <b>A</b>ll <b>L</b>eptons.
 *
 * \section intro_sec Introduction
 *
 * This is the introduction.
 *
 * \section install_sec Installation
 *
 * \subsection requirements
 *
 * \section HowTo
 *
 */



/**
 * main mmc class - must be constructed before anything else
 */

class Propagate: public PhysicsModel
{

protected:

    PROPOSALParticle *particle_;
    Medium *medium_;
    CrossSections *cros;
    Energy2Loss *E2Loss_;
    Output *o;
    double rho;		//!< multiplicative medium density correction factor

    StandardNormal *StandardN;
    std::vector<Integral *> integral_;

    bool pint;
    double decayS, ionizS, bremsS, epairS, photoS, totalS;

    double ec, tc;

    /*!
     * \brief indicates if the interpolated function is increasing or decreasing.
     *
     * True  =   Interpolated function is increasing. \n
     * False =   Interpolated function is decreasing. BigLow is set to f(x_min)
     *           and the interpolation points that are saved are BigLow-f(x).
    */
    bool up;

    double bigLow[2];               //!< See the describtion of the variable up.
    double storeDif[2];             //!< Stores the interpolated values for farther calculations

    bool df;                        //!< Differential. Chooses if function is calculated is return differential. Preset as false in constructor
    Interpolate *interpolateJ_;     //!< Interpolate object of the Integral of the function
    Interpolate *interpolateJdf_;   //!< Interpolate object of the function


public:

//----------------------------------------------------------------------------//

    bool jt;            //!< True = Interpolation is activated. Preset as false in constructor
    bool recc;          //!< Set as false in constructor
    bool exactTime;     //!< True = Interpolate the time difference of the advanced particle. Set as false in constructor
    bool contiCorr;     //!< True = Continious Correction is activated. Preset as false in constructor
    bool molieScat;     //!< True = Molie Scattering is activated. Preset as false in constructor
    bool sdec;          //!< True = particle decayed during propagation
    static int g;		//!< Order of interpolation
    bool dw;            //!< Do weigthing? Set as false in constructor
    double rw;	      	//!< Re-weighting order. Set to 0 in constructor
    double hw;	      	//!< Distance at which re-weighting starts. Set to 0 in constructor

//----------------------------------------------------------------------------//
    // Constructors
    /**
     * empty default constructor.
     */
    Propagate();


//----------------------------------------------------------------------------//
    /**
     * initialize all classes necessary for propagation of a muon.
     *
     * \param   w       medium to propagate through
     * \param   ecut    stochastic calculation when energy > ecut
     * \param   vcut    stochastic calculation when e-loss > energy*vcut
     */
    Propagate(std::string w, double ecut, double vcut);

//----------------------------------------------------------------------------//
    /**
     * initialize all classes necessary for propagation of a muon or tau.
     *
     * \param   w       medium to propagate through
     * \param   ecut    stochastic calculation when energy > ecut
     * \param   vcut    stochastic calculation when e-loss > energy*vcut
     * \param   type    particle type
     * \return
     */

    Propagate(std::string w, double ecut, double vcut, std::string type);

//----------------------------------------------------------------------------//
    /**
     * initialize all classes necessary for propagation of a muon or tau.
     *
     * \param   w       medium to propagate through
     * \param   ecut    stochastic calculation when energy > ecut
     * \param   vcut    stochastic calculation when e-loss > energy*vcut
     * \param   type    particle type
     * \param   rho     multiplicative medium density correction factor
     */

    Propagate(std::string w, double ecut, double vcut, std::string type, double rho);

//----------------------------------------------------------------------------//

    /**
     * class initializer
     *
     * \param   w       medium to propagate through
     * \param   ecut    stochastic calculation when energy > ecut
     * \param   vcut    stochastic calculation when e-loss > energy*vcut
     * \param   type    particle type
     * \param   rho     multiplicative medium density correction factor
     */
    void init(std::string w, double ecut, double vcut, std::string type, double rho);

//----------------------------------------------------------------------------//
    /**
     * Propagates the particle of initial energy e to the distance r.
     * Returns the final energy if the
     * particle has survived or the track length to the
     * point of disappearance with a minus sign otherwise.
     *
     *  \image html PAP-PROPOSAL.jpg
     *  \image latex PAP-PROPOSAL.pdf "PAP of the propagateTo routine" width=10cm
     *
     *  \param  r   maximum track length
     *  \param  e   initial energy
     *  \return energy at distance r OR -(track length)
     */

    double propagateTo(double r, double e);

//----------------------------------------------------------------------------//
    /**
     * Propagates the particle of initial energy e to the distance r.
     * Returns the final energy if the
     * particle has survived or the track length to
     * the point of disappearance with a minus sign otherwise.
     * Also calculates particle energy at point rc.
     * Call getPropEc() to get this energy.
     *
     *  \param  r   maximum track length
     *  \param  e   initial energy
     *  \param  rc  track length at which the energy of the particle
     *              is additional calculated
     *  \return energy at distance r OR -(track length)
     */

    double propagateTo(double r, double e, double rc);

//----------------------------------------------------------------------------//
    /**
     * returns the particle energy at rc if the particle has survived or the
     * distance from the point of decay to rc with a minus sign otherwise.
     *
     *  \return energy at distance rc OR -(|distance of decay - rc|)
     */

    double getPropEc();

//----------------------------------------------------------------------------//
    /**
     * returns the particle time at rc if the particle
     * has survived or at the point of decay otherwise.
     *
     *  \return particle time at distance rc OR particle time at the decay
     */

    double getPropTc();

//----------------------------------------------------------------------------//
    /**
     * function for energy range calculation - interface to Integral
     *
     *  \param  E   particle energy
     *  \return Returns the probability [1/MeV]
     */

    double function(double E);

//----------------------------------------------------------------------------//
    /**
     * call this routine to enable interpolations.
     * To enable everything, set w="all"
     *
     *  \param  w   defines which integrals should be interpolated
     */

    void interpolate(std::string w);

//----------------------------------------------------------------------------//

    /**
     * call this routine to enable interpolations and
     * their save and reread. To enable everything, set w="all"
     *
     *  \param  w           defines which integrals should be interpolated
     *  \param  filename    file in which the interpolation tables will be saved
     *                      or read from. "" will create a default file name
     */

    void interpolate(std::string w, std::string filename);

//----------------------------------------------------------------------------//

    /**
     * Calculates the value of the tracking integral
     *
     *  \param  ei      initial energy
     *  \param  rnd     random number which is used for the calculation
     *  \param  pint    particle interaction? (false = decay)
     *  \return value of the tracking integral [ 1 ]
     */

    double getpr(double ei, double rnd, bool pint);

//----------------------------------------------------------------------------//

    /**
     * final energy, corresponding to the given value rnd of the
     * tracking integral
     *
     *  \param  ei      initial energy
     *  \param  rnd     random number which is used for the calculation
     *  \param  pint    particle interaction? (false = decay)
     *  \return final energy due to continous energy losses [MeV]
     */

    double getef(double ei, double rnd, bool pint);

//----------------------------------------------------------------------------//
    /**
     * 1d parametrization - interface to Interpolate
     *
     *  \param  e   particle energy
     *  \return integrated or differential function value
     */

    double functionInt(double e);
//----------------------------------------------------------------------------//

    /**
     * Set cross section parametrizations
     */

    void setCSform(int bspar, int pncrs, int pncbb, int pncsh, std::string muta);
//----------------------------------------------------------------------------//

    /**
     * Set cross sections to highest parametrizations
     */

    void setCShigh();
//----------------------------------------------------------------------------//

    /**
     * Set cross sections to lowest parametrizations
     */

    void setCSlow();
//----------------------------------------------------------------------------//
    // Getter

    /**
     * Getter for the exactTime.
     *
     *  \return exactTime is calculated?
     */
    bool get_exactTime()
    {
        return exactTime;
    }

    /**
     * Getter for molieScat.
     *
     *  \return molie scattering is activated?
     */
    bool get_molieScat()
    {
        return molieScat;
    }

    /**
     * Getter for rho.
     *
     *  \return multiplicative factor for the medium density
     */
    double get_rho()
    {
        return rho;
    }

    StandardNormal *get_Standard();

    CrossSections *get_cros();

    Medium *get_Medium();

    PROPOSALParticle *get_particle()
    {
        return particle_;
    }

    Output *get_output()
    {
        return o;
    }
//----------------------------------------------------------------------------//
    //Setter

    void set_exactTime(bool eTime);

    void set_molieScat(bool mScat);

    void set_contiCorr(bool cCorr);

    void set_pint(bool newPint);

    void set_rho(double newRho);
    
    virtual void SetRandomNumberGenerator(boost::function<double ()> &);

};


#endif /* PROPAGATE_H_ */
