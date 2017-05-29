/*! \file   Scattering.h
*   \brief  Header file for the Scattering bug routines.
*
*   This version has a major bug and produces too small scattering angles.
*
*   \date   2013.08.19
*   \author Tomasz Fuchs
*/


#ifndef SCATTERING_H
#define SCATTERING_H

#include <vector>
#include <string>

#include "PROPOSAL/PROPOSALParticle.h"
#include "PROPOSAL/CrossSections.h"
#include "PROPOSAL/Interpolant.h"
#include "PROPOSAL/Integral.h"
#include "PROPOSAL/MathModel.h"

namespace PROPOSAL{

/**
  * \brief This class provides the scattering routine provided by moliere.
  *
  * More precise scattering angles will be added soon.
  */


class Scattering : public MathModel
{


private:

    double x0_;

    bool do_interpolation_;
    int order_of_interpolation_;

    Integral* integral_;
    Interpolant* interpolant_;
    Interpolant* interpolant_diff_;

    PROPOSALParticle* particle_;
    PROPOSALParticle* backup_particle_;
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

    void SetParticle(PROPOSALParticle* particle);
    void SetCrosssections(std::vector<CrossSections*> crosssections);
//----------------------------------------------------------------------------//
    // Getter

    PROPOSALParticle* GetParticle()
    {
        return particle_;
    }

    double GetX0()
    {
        return x0_;
    }

//----------------------------------------------------------------------------//
    // destructors
    ~Scattering() {}

    PROPOSALParticle *GetBackup_particle() const;
    void SetBackup_particle(PROPOSALParticle *backup_particle);
    void RestoreBackup_particle();

};
};


#endif //SCATTERING_H
