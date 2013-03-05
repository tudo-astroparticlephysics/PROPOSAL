/*! \file   Mpairproduction.h
*   \brief  Header file for the pairproduction routines.
*
*   For more details see the class documentation.
*
*   \date   15.03.2011
*   \author Martin Schmitz
*/

#ifndef MPAIRPRODUCTION_H
#define MPAIRPRODUCTION_H


#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Interpolate.h"



/*! \class Mpairproduction Mpairproduction.h "Mpairproduction.h"
 *  \brief class contains functions necessary for calculation of
 *  \mu+/\mu- pair production losses.
 *  All parametrizations are from Kelner, S., Kokoulin, R., & Petrukhin, A.,
 *  Direct Production of Muon Pairs by High-Energy Muons,Phys.
 *  of Atomic Nuclei, Vol. 63, No 9 (2000) 1603.
 *  All variable names are consinstent with this paper.
 */

class CrossSections;
class MpairContinuous;
class MpairStochastic;

class Mpairproduction: public CrossSections
{

protected:

    Integral*                   integral_;        ///< Object for integration
    std::vector<Interpolate*>   interpolateJ_;    ///< Object for Interpolation

    Mpairproduction *Mpair_;                    ///< Pointer to remember everything
    MpairContinuous *continuous_;               ///< Pointer to remember everything
    MpairStochastic *stochastic_;               ///< Pointer to remember everything

    double v;                                   ///< energy fraction transferred to the particles of a pair
    double vMin_;                               ///< integrationlimit for v integration
    double vUp_;                                ///< integrationlimit for v integration
    double vMax_;                               ///< integrationlimit for v integration
    double rhomax;                              ///< integrationlimit for \rho integration

    int NumberOfComponent;                      ///< Which Component of media is used for mpairproduction
    int form_;                                  ///< parametrizationswitcher
    int cm_;                                    ///< turn mpair off with cm = 0;


//----------------------------------------------------------------------------//

    /**
    * Sets the paramter for integration/interpolation
    *
    * \param    i   crossections component
    */
    void setEnergy(int i);
//----------------------------------------------------------------------------//

    /**
    * function which will be integrated.
    *
    * \param    x   rho???
    * \return   ???
    */
    double function(double x);
//----------------------------------------------------------------------------//

    /**
    * function which will be interpolated
    *
    * \param    e   particle energy
    * \param    v   relative energy loss
    * \return   ???
    */
    using CrossSections::functionInt;
    double functionInt(double e, double v);
//----------------------------------------------------------------------------//

    /**
     * This is what U is equal to. See eq. (22) in
     * Kelner, S., Kokoulin, R., & Petrukhin, A.,
     * Direct Production of Muon Pairs by High-Energy Muons,Phys.
     * of Atomic Nuclei, Vol. 63, No 9 (2000) 1603 for more information
     *
     * \param   rho
     * \return  U
     */
    double calculateU(double rho);
//----------------------------------------------------------------------------//


public:

    /**
    * usual Constructor
    *
    * \param    *cros   crossection which will be used for energy losses
    */
    Mpairproduction(CrossSections *cros);

//----------------------------------------------------------------------------//

    /**
     * Copy Constructor from other Mpairproduction class.
     *
     * \param    *cros   source class which will be copied.
     */
    Mpairproduction(Mpairproduction *cros_);

//----------------------------------------------------------------------------//
    /**
     * this is what \simga dv is equal to
     *
     * \param   v   relative energy loss
     * \param   i   component number
     * \return  ???
     */

    double mpair(double v, int i);
//----------------------------------------------------------------------------//

    /**
     * \brief activates interpolation
     */
    void activate_interpolation();

//----------------------------------------------------------------------------//
    // Getter

    MpairContinuous* get_Continuous()
    {
        return continuous_;
    }

    int get_cm()
    {
        return cm_;
    }

    MpairStochastic* get_Stochastic()
    {
        return stochastic_;
    }

//----------------------------------------------------------------------------//
    // Setter

    void set_cm(int newcm)
    {
        cm_ =   newcm;
    }


};

#endif // MPAIRPRODUCTION_H
