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

ScatteringNoScattering::ScatteringNoScattering(const ParticleDef& particle_def, std::shared_ptr<const Medium> medium)
    : Scattering(particle_def)
    , medium_(medium)
{
}

ScatteringNoScattering::ScatteringNoScattering(const ScatteringNoScattering& scattering)
    : Scattering(scattering)
    , medium_(scattering.medium_)
{
}

ScatteringNoScattering::ScatteringNoScattering(const ParticleDef& particle_def, const ScatteringNoScattering& scattering)
    : Scattering(particle_def)
    , medium_(scattering.medium_)
{
}

ScatteringNoScattering::~ScatteringNoScattering()
{ }

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

void ScatteringNoScattering::print(std::ostream& os) const
{
    os << "Medium:\n" << *medium_ << '\n';
}

//----------------------------------------------------------------------------//
Scattering::RandomAngles ScatteringNoScattering::CalculateRandomAngle(double dr,
                                                                      double ei,
                                                                      double ef,
                                                                      const Vector3D& pos,
                                                                      double rnd1,
                                                                      double rnd2,
                                                                      double rnd3,
                                                                      double rnd4)
{
    (void)ei;
    (void)ef;
    (void)dr;
    (void)pos;
    (void)rnd1;
    (void)rnd2;
    (void)rnd3;
    (void)rnd4;

    Scattering::RandomAngles random_angles;

    random_angles.sx = 0;
    random_angles.tx = 0;

    random_angles.sy = 0;
    random_angles.ty = 0;

    return random_angles;
}
