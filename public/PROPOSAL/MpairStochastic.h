/*! \file   MpairStochastic.h
*   \brief  Header file for the stochastic pairproduction routines.
*
*   \author Jan-Hendrik Koehne
*/
#ifndef MPAIRSTOCHASTIC_H
#define MPAIRSTOCHASTIC_H

#include "PROPOSAL/Mpairproduction.h"

class MpairStochastic : public Mpairproduction {

protected:

    std::vector<Integral*>      integral_;
    std::vector<double>         H_;
    std::vector<Interpolate*>   interpolateJ_;
    std::vector<Interpolate*>   interpolateJo_;
    double sum;
    double rnd;

//----------------------------------------------------------------------------//
    /*!
     * \brief   ??????
     *
     * \param   v   relative energy loss
     * \return  ???
     */
    double functionInt(double v);

//----------------------------------------------------------------------------//
    /*!
     * \brief   ??????
     *
     * \param   e   particle energy
     * \param   v   relative energy loss
     * \return  ???
     */
    double functionInt(double e, double v);

//----------------------------------------------------------------------------//
    /*!
     * \brief   ??????
     *
     * \param   v   relative energy loss
     * \return  ???
     */
    double function(double v);

//----------------------------------------------------------------------------//

public:
    /*!
     * \brief default constructor
     *
     * \param   *MPAIR  crossection for pairproduction
     */
    MpairStochastic(Mpairproduction *MPAIR);

//----------------------------------------------------------------------------//

    /*!
     * \brief   ??????
     *
     * \return  ???
     */
    double dNdx();
//----------------------------------------------------------------------------//

    /*!
     * \brief   ??????
     *
     * \param   rnd     random number
     * \return  ???
     */
    double dNdx(double rnd);

//----------------------------------------------------------------------------//
    /*!
     * \brief   ??????
     *
     * \param   rnd     random number
     * \return  ???
     */
    double e(double rnd);

//----------------------------------------------------------------------------//
    /**
     * \brief activates interpolation
     */
    void activate_interpolation();


};

#endif // MPAIRSTOCHASTIC_H
