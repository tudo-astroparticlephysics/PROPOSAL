#pragma once

#include <boost/bimap.hpp>

#include "PROPOSAL/PROPOSAL.h"
#include "dataclasses/physics/I3Particle.h"

namespace I3PROPOSALParticleConverter {

I3Particle::ParticleType GenerateI3Type(const PROPOSAL::DynamicData& secondary);

PROPOSAL::ParticleDef GeneratePROPOSALType(const I3Particle::ParticleType& ptype_I3);

} // namespace I3PROPOSALParticleConverter
