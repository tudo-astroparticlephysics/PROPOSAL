/*! \file   PhysicsModel.h
*   \brief  Header file for the physics model routines.
*
*   For more details see the class documentation.
*
*   \date   21.06.2010
*   \author Jan-Hendrik Koehne
*/

#ifndef PhysicsModel_H
#define PhysicsModel_H

#include <cmath>
#include "PROPOSAL/FunctionInt.h"
#include "PROPOSAL/FunctionOfx.h"
#include "PROPOSAL/FunctionInt2.h"
#include "PROPOSAL/MathModel.h"



//#ifndef
//    #define
//#endif

/*! \class PhysicsModel PhysicsModel.h "PhysicsModel.h"
 *  \brief This is an entry class of the physical model we are building.
 */

class PhysicsModel: public FunctionInt, public FunctionOfx, public FunctionInt2, public MathModel
{


protected:

    static double elow_;			// bounds of parameterizations
    static double nlow_;				// bounds of parameterizations
    static double ebig_;			// bounds of parameterizations




public:

//----------------------------------------------------------------------------//
    // constructors

    /**
     * \brief  Create a PhysicsModel.
     */
    PhysicsModel();


//----------------------------------------------------------------------------//
    // Memberfunctions

    /**
     * \brief  for integration - interface to Integral - to be overwritten by subclasses
     *
     * \param   e   dummy parameter
     * \return  0
     */

    double function(double e)
    {
        return 0;
    }

//----------------------------------------------------------------------------//

    /**
     * \brief  1d parametrization - interface to Interpolate - to be overwritten by subclasses
     *
     * \param   e   dummy parameter
     * \return  0
     */
    double functionInt(double e)
    {
        return 0;
    }
//----------------------------------------------------------------------------//

    /**
     * \brief  2d parametrization - interface to Interpolate - to be overwritten by subclasses
     *
     * \param   e   dummy parameter
     * \param   v   dummy parameter
     * \return  0
     */
    double functionInt(double e, double v)
    {
        return 0;
    }
//----------------------------------------------------------------------------//
    // Getter

    static double get_ebig()
    {
        return ebig_;
    }

    static double get_elow()
    {
        return elow_;
    }
//----------------------------------------------------------------------------//
    // Setter

    static void set_ebig( double newebig )
    {
        ebig_ = newebig;
    }

    static void set_elow( double newelow )
    {
        elow_ = newelow;
    }

//----------------------------------------------------------------------------//
    // destructors

    /**
     * \brief Crush this PhysicsModel.
     */
    virtual ~PhysicsModel();


};



#endif //PhysicsModel_H
