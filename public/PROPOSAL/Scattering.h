/*! \file   Scattering.h
*   \brief  Header file for the Scattering routines.
*
*   For more details see the class documentation.
*
*   \date   2013.06.13
*   \author Tomasz Fuchs
*/


#ifndef SCATTERING_H
#define SCATTERING_H
#include "vector"
#include <string>
#include "PROPOSAL/ProcessCollection.h"
#include "PROPOSAL/Particle.h"
#include "PROPOSAL/Interpolant.h"
#include "PROPOSAL/Integral.h"

/**
  * \brief This class provides the scattering routine provided by moliere.
  *
  * More precise scattering angles will be added soon.
  */


class Scattering
{


private:

    double x0_; //Radiation Length
    bool do_interpolation_;
    int order_of_interpolation_;

    Integral* integral_;
    Interpolant* interpolant_;
    Interpolant* interpolant_diff_;

    Particle* particle_;
    std::vector<CrossSections*> crosssections_;
//----------------------------------------------------------------------------//

public:

    /**
     * \brief Default Constructor
     *
     * Constructor which sets "default" settings.
     */
    Scattering();

//----------------------------------------------------------------------------//

    Scattering(std::vector<CrossSections*> crosssections);
//----------------------------------------------------------------------------//

    Scattering(const Scattering&);
    Scattering& operator=(const Scattering&);
    bool operator==(const Scattering &scattering) const;
    bool operator!=(const Scattering &scattering) const;
//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
    // Memberfunctions

    double  CalculateTheta0(double dr, double ei, double ef);
    double  FunctionToIntegral(double energy);
    double  FunctionToBuildInterpolant(double energy);
    void    EnableInterpolation();
    void    DisableInterpolation();
//----------------------------------------------------------------------------//

    void swap(Scattering &scattering);

//----------------------------------------------------------------------------//

    //Setter

    void SetParticle(Particle* particle);
    void    SetCrossSections(std::vector<CrossSections*> crosssections);
//----------------------------------------------------------------------------//
    // Getter


//----------------------------------------------------------------------------//
    // destructors
    ~Scattering() {}


};



#endif //SCATTERING_H
