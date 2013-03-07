/*! \file   MpairContinuous.h
*   \brief  Header file for the continuous pairproduction routines.
*
*   \author Jan-Hendrik Koehne
*/

#ifndef MPAIRCONTINUOUS_H
#define MPAIRCONTINUOUS_H
#include"PROPOSAL/Mpairproduction.h"

class MpairContinuous : public Mpairproduction{


protected:

    Interpolate *interpolateJ_;
//----------------------------------------------------------------------------//

    /*!
     * \brief   ??????
     *
     * \param   v   relative energy loss
     * \return  ???
     */
    double function(double v);
//----------------------------------------------------------------------------//

    /*!
     * 1d parametrization - interface to Interpolate
     *
     * \param   e   energy
     * \return  dE/dx for particle energy e
     */
    using Mpairproduction::functionInt;
    double functionInt(double e);

//----------------------------------------------------------------------------//
public:

    /*!
     * \brief default constructor
     *
     * \param   *MPAIR  crossection for pairproduction
     */

    MpairContinuous(Mpairproduction *MPAIR);
//----------------------------------------------------------------------------//
    /*!
     * \brief This calculates the contribution to dEdx() up to v_max.
     *
     * \return dE/dx [MeV/cm]
     */
    double dEdx();

//----------------------------------------------------------------------------//
    /**
     * \brief activates interpolation
     */
    void activate_interpolation();



};

#endif // MPAIRCONTINUOUS_H
