/*! \file   PhotoContinuous.h
*   \brief  Header file for the continuous photonuclear routines.
*
*   For more details see the class documentation.
*
*   \date   28.06.2010
*   \author Jan-Hendrik Koehne
*/

#ifndef PHOTOCONTINUOUS_H_
#define PHOTOCONTINUOUS_H_

#include "PROPOSAL/Photonuclear.h"
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Interpolate.h"



/*! \class PhotoContinuous PhotoContinuous.h "PhotoContinuous.h"
   \brief class contains functions necessary for calculation of photonuclear losses
 */




class PhotoContinuous: public Photonuclear{

protected:

    Integral *integral_;
    bool jt_;

//----------------------------------------------------------------------------//


public:

    Interpolate *interpolateJ_;
//----------------------------------------------------------------------------//
    //Constructor

    /*!
     * \brief initializes an integral class separate from that in Propagate
     *
     * \param   *cros   photonuclear crossection
     */

    PhotoContinuous(Photonuclear *cros);

//----------------------------------------------------------------------------//
    //Memberfunctions

    /*!
     * \brief photonuclear energy losses - interface to Integral
     *
     * \param   v   relative energy loss
     * \return  \f[v\cdot photoN\f]
    */

    double function(double v);

//----------------------------------------------------------------------------//

    /*!
     * \brief Contribution of photonuclear interactions to -dE/dx
     *
     * \f[\frac{dE}{dx} = c_p E_p \sum_{numC}
     * \int_{v_{Min}}^{v_{Up}} v \frac{d\sigma}{dv} dv\f]
     *
     * \return dE/dx [MeV/cm]
     */

    double dEdx();

//----------------------------------------------------------------------------//
    /*!
     * 1d parametrization - interface to Interpolate
     *
     * \param   e   sets particle to energy e
     * \return  dE/dx
     */
    using Photonuclear::functionInt;
    double functionInt(double e);

//----------------------------------------------------------------------------//
    void set_jt(bool newJT)
    {
        jt_ =   newJT;
    }

//----------------------------------------------------------------------------//
    // Destructor

    ~PhotoContinuous(){}

};


#endif /* PHOTOCONTINUOUS_H_ */
