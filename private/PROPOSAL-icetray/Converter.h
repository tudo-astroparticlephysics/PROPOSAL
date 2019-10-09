#pragma once

#include <boost/bimap.hpp>

#include "PROPOSAL/PROPOSAL.h"
#include "dataclasses/physics/I3Particle.h"

class I3PROPOSALParticleConverter {
public:
    I3PROPOSALParticleConverter();

    /**
     * Get the internal MMC name associated with a particle type
     */
    I3Particle::ParticleType GenerateI3Type(const PROPOSAL::DynamicData& secondary) const;

    const PROPOSAL::ParticleDef& GeneratePROPOSALType(const I3Particle::ParticleType& ptype_I3) const;

    PROPOSAL::Particle GeneratePROPOSALParticle(const I3Particle& i3_particle) const;

    I3Particle GenerateI3Particle(const PROPOSAL::DynamicData& proposal_particle) const;
private:
    std::map<I3Particle::ParticleType, const PROPOSAL::ParticleDef&> i3_to_proposal_;
    PROPOSAL::ParticleDef null_definition_;

};
