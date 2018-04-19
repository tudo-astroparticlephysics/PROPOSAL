/*! \file   ScatteringNoScattering.cxx
 *   \brief  Source file for the ScatteringNoScattering routines.
 *
 *   For more details see the class documentation.
 *
 *   \date
 *   \author
 */

#include "PROPOSAL/scattering/ScatteringNoScattering.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/Particle.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

ScatteringNoScattering::ScatteringNoScattering(Particle& particle, const Medium& medium)
    : Scattering(particle)
    , medium_(medium.clone())
{
}

ScatteringNoScattering::ScatteringNoScattering(const ScatteringNoScattering& scattering)
    : Scattering(scattering)
    , medium_(scattering.medium_->clone())
{
}

ScatteringNoScattering::ScatteringNoScattering(Particle& particle, const ScatteringNoScattering& scattering)
    : Scattering(particle)
    , medium_(scattering.medium_->clone())
{
}

ScatteringNoScattering::~ScatteringNoScattering()
{
    delete medium_;
}

bool ScatteringNoScattering::compare(const Scattering& scattering) const
{
    const ScatteringNoScattering* scatteringNoScattering = dynamic_cast<const ScatteringNoScattering*>(&scattering);

    if (!scatteringNoScattering)
        return false;
    else if (*medium_ != *scatteringNoScattering->medium_)
        return false;
    else
        return true;
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
