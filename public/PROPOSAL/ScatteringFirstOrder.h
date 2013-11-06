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
#include "PROPOSAL/MathModel.h"
#include "PROPOSAL/Medium.h"

/**
  * \brief This class provides the scattering routine provided by moliere.
  *
  * More precise scattering angles will be added soon.
  */


class ScatteringFirstOrder : public MathModel
{

//----------------------------------------------------------------------------//

public:

    /**
     * \brief Default Constructor
     *
     * Constructor which sets "default" settings.
     */
    ScatteringFirstOrder();


//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
    // Memberfunctions

    double  CalculateTheta0(double dr, Particle* part, Medium* med);
    void    Scatter(double dr, Particle* part, Medium* med);


//----------------------------------------------------------------------------//
    // destructors
    ~ScatteringFirstOrder() {}


};



#endif //SCATTERING_FIRSTORDER_H
