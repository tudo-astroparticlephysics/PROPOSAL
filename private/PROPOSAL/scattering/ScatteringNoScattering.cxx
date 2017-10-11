/*! \file   ScatteringNoScattering.cxx
*   \brief  Source file for the ScatteringNoScattering routines.
*
*   For more details see the class documentation.
*
*   \date   
*   \author 
*/

#include "PROPOSAL/particle/PROPOSALParticle.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/scattering/ScatteringNoScattering.h"


using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //


ScatteringNoScattering::ScatteringNoScattering(PROPOSALParticle& particle, const Medium& medium)
    : Scattering(particle)
    , medium_(medium.clone())
{
}

ScatteringNoScattering::ScatteringNoScattering(const ScatteringNoScattering& scattering)
    : Scattering(scattering)
    , medium_(scattering.medium_->clone())
{
}

ScatteringNoScattering::ScatteringNoScattering(PROPOSALParticle& particle, const ScatteringNoScattering& scattering)
    : Scattering(particle)
    , medium_(scattering.medium_->clone())
{
}

ScatteringNoScattering::~ScatteringNoScattering()
{
    delete medium_;
}


//----------------------------------------------------------------------------//

Scattering::RandomAngles ScatteringNoScattering::CalculateRandomAngle(double dr, double ei, double ef)
{
    (void)ei;
    (void)ef;
    (void)dr;

    Scattering::RandomAngles random_angles;

    random_angles.sx = 0;
    random_angles.tx = 0;

    random_angles.sy = 0;
    random_angles.ty = 0;

    return random_angles;
}
