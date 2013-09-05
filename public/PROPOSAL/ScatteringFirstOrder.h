/*! \file   Scattering.h
*   \brief  Header file for the Scattering routines.
*
*   For more details see the class documentation.
*
*   \date   2013.06.13
*   \author Tomasz Fuchs
*/


#ifndef SCATTERING_FIRSTORDER_H
#define SCATTERING_FIRSTORDER_H
#include "vector"
#include <string>
#include "PROPOSAL/Particle.h"
#include "PROPOSAL/StandardNormal.h"
#include "PROPOSAL/Medium.h"

/**
  * \brief This class provides the scattering routine provided by moliere.
  *
  * More precise scattering angles will be added soon.
  */


class ScatteringFirstOrder : public MathModel
{


private:

    StandardNormal* standard_normal_;
//----------------------------------------------------------------------------//

public:

    /**
     * \brief Default Constructor
     *
     * Constructor which sets "default" settings.
     */
    ScatteringFirstOrder();

    ScatteringFirstOrder(StandardNormal* standard_normal);

//----------------------------------------------------------------------------//

    ScatteringFirstOrder(const ScatteringFirstOrder&);
    ScatteringFirstOrder& operator=(const ScatteringFirstOrder&);
    bool operator==(const ScatteringFirstOrder &scattering) const;
    bool operator!=(const ScatteringFirstOrder &scattering) const;
//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
    // Memberfunctions

    double  CalculateTheta0(double dr, Particle* part, Medium* med);
    void    Scatter(double dr, Particle* part, Medium* med);

//----------------------------------------------------------------------------//

    void swap(ScatteringFirstOrder &scattering);

//----------------------------------------------------------------------------//

    StandardNormal* GetStandardNormal()
    {
        return standard_normal_;
    }

//----------------------------------------------------------------------------//
    // destructors
    ~ScatteringFirstOrder() {}


};



#endif //SCATTERING_FIRSTORDER_H
