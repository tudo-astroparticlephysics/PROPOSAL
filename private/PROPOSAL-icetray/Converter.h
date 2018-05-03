#pragma once

#include <boost/bimap.hpp>

#include "PROPOSAL/PROPOSAL.h"
#include "dataclasses/physics/I3Particle.h"

namespace I3PROPOSALParticleConverter {

/**
 * Get the internal MMC name associated with a particle type
 */
I3Particle::ParticleType GenerateI3Type(const PROPOSAL::DynamicData& secondary);

PROPOSAL::ParticleDef GeneratePROPOSALType(const I3Particle::ParticleType& ptype_I3);

PROPOSAL::Particle GeneratePROPOSALParticle(const I3Particle& i3_particle);

I3Particle GenerateI3Particle(const PROPOSAL::DynamicData& proposal_particle);

} // namespace I3PROPOSALParticleConverter
