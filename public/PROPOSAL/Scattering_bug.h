/*! \file   Scattering_bug.h
*   \brief  Header file for the Scattering bug routines.
*
*   This version has a major bug and produces too small scattering angles.
*
*   \date   2013.08.19
*   \author Tomasz Fuchs
*/


#ifndef SCATTERING_BUG_H
#define SCATTERING_BUG_H
#include "vector"
#include <string>
#include "PROPOSAL/Particle.h"
#include "PROPOSAL/CrossSections.h"
#include "PROPOSAL/StandardNormal.h"
#include "PROPOSAL/Interpolant.h"
#include "PROPOSAL/Integral.h"

/**
  * \brief This class provides the scattering routine provided by moliere.
  *
  * More precise scattering angles will be added soon.
  */


class Scattering : public MathModel
{


private:

    bool do_interpolation_;
    int order_of_interpolation_;

    Integral* integral_;
    Interpolant* interpolant_;
    Interpolant* interpolant_diff_;

    Particle* particle_;
    std::vector<CrossSections*> crosssections_;
    StandardNormal* standard_normal_;
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

    long    double  CalculateTheta0(double dr, double ei, double ef);
    void            Scatter(double dr, double ei, double ef);
    double          FunctionToIntegral(double energy);
    double          FunctionToBuildInterpolant(double energy);
    void            EnableInterpolation(std::string path = "");
    void            DisableInterpolation();
//----------------------------------------------------------------------------//

    void swap(Scattering &scattering);

//----------------------------------------------------------------------------//

    //Setter

    void SetParticle(Particle* particle);
    void SetCrosssections(std::vector<CrossSections*> crosssections);
//----------------------------------------------------------------------------//
    // Getter

    Particle* GetParticle()
    {
        return particle_;
    }

//----------------------------------------------------------------------------//
    // destructors
    ~Scattering() {}


};



#endif //SCATTERING_BUG_H
