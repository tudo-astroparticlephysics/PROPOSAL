/*! \file   BremsContinuous.h
*   \brief  Header file for definition of the continuous-Bremsstrahlung class object.
*
*   For more details see the class documentation.
*
*   \date   21.06.2010
*   \author Jan-Hendrik Koehne
*/

#ifndef BremsContinuous_H
#define BremsContinuous_H

#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/Interpolate.h"
#include <cmath>

/*! \class BremsContinuous BremsContinuous.h "BremsContinuous.h"
    \brief class contains functions necessary
           for calculation of bremsstrahlung losses
 */

class BremsContinuous: public Bremsstrahlung
{

public:

    Interpolate *interpolateJ_;

    // constructors

//----------------------------------------------------------------------------//

    /*!
    Create a BremsContinuous.
    */

    BremsContinuous();

//----------------------------------------------------------------------------//

    /*!
    initializes an integral class separate from that in Propagate
    */

    BremsContinuous(Bremsstrahlung *cros);

//----------------------------------------------------------------------------//

    // Memberfunctions

    /*!
    bremsstrahlung energy losses - interface to Integral;
    \return \f[ v \cdot \sigma\; dv \f]
    */

    double function(double v);

//----------------------------------------------------------------------------//

    /*!
    contribution of bremsstrahlung to -dE/dx;
    \f[ \frac{dE}{dx} =
        c_b E_p \sum_{numC} \int_{v_{Min}}^{v_{Up}} v \sigma\; dv\f]
    \return dEdx [E/cm]
    */

    double dEdx();

//----------------------------------------------------------------------------//

    /*!
    1d parametrization - interface to Interpolate
    */

    double functionInt(double e);

//----------------------------------------------------------------------------//

    // Getter

    Integral* get_integral() const
    {
        return integral_;
    }

    Interpolate* get_interpolateJ()
    {
        return interpolateJ_;
    }

//----------------------------------------------------------------------------//

    // Setter

    void set_integral(Integral* integral);
    void set_interpolateJ(Interpolate* interpolateJ);
    void set_jt(bool newjt);

//----------------------------------------------------------------------------//

    // destructors

    ///@brief Crush this BremsContinuous.

    virtual ~BremsContinuous();

};


#endif //BremsContinuous_H
