
#include <boost/assign.hpp>

#include "Converter.h"

namespace I3PROPOSALParticleConverter {

typedef boost::bimap<I3Particle::ParticleType, std::string> bimap_ParticleType;
static const bimap_ParticleType I3_PROPOSAL_ParticleType_bimap =
    boost::assign::list_of<bimap_ParticleType::relation>(I3Particle::MuMinus, "MuMinus")(I3Particle::MuPlus, "MuPlus")(
        I3Particle::TauMinus,
        "TauMinus")(I3Particle::TauPlus, "TauPlus")(I3Particle::EMinus, "EMinus")(I3Particle::EPlus, "EPlus")(
        I3Particle::NuMu,
        "NuMu")(I3Particle::NuMuBar, "NuMuBar")(I3Particle::NuE, "NuE")(I3Particle::NuEBar, "NuEBar")(
        I3Particle::NuTau,
        "NuTau")(I3Particle::NuTauBar, "NuTauBar")(I3Particle::Brems, "Brems")(I3Particle::DeltaE, "DeltaE")(
        I3Particle::PairProd,
        "EPair")(I3Particle::NuclInt, "NuclInt")(I3Particle::MuPair, "MuPair")(I3Particle::Hadrons, "Hadrons")(
        I3Particle::Monopole,
        "Monopole")(I3Particle::STauMinus, "STauMinus")(I3Particle::STauPlus, "STauPlus")(I3Particle::Gamma, "Gamma")(
        I3Particle::Pi0,
        "Pi0")(I3Particle::PiPlus, "PiPlus")(I3Particle::PiMinus, "PiMinus")(I3Particle::K0_Short,
                                                                             "K0") // TODO(mario):  Fri 2017/11/17
    (I3Particle::KPlus, "KPlus")(I3Particle::KMinus, "KMinus")(I3Particle::PPlus, "PPlus")(I3Particle::PMinus,
                                                                                           "PMinus");

} // namespace I3PROPOSALParticleConverter

// ------------------------------------------------------------------------- //
I3Particle::ParticleType I3PROPOSALParticleConverter::GenerateI3Type(const PROPOSAL::DynamicData& secondary)
{
    using namespace PROPOSAL;

    DynamicData::Type type = secondary.GetTypeId();

    switch (type)
    {
        case DynamicData::Particle:
        {
            const PROPOSAL::Particle& particle = static_cast<const PROPOSAL::Particle&>(secondary);
            ParticleDef particle_def           = particle.GetParticleDef();

            I3Particle::ParticleType ptype_I3;

            bimap_ParticleType::right_const_iterator proposal_iterator =
                I3_PROPOSAL_ParticleType_bimap.right.find(particle_def.name);
            if (proposal_iterator == I3_PROPOSAL_ParticleType_bimap.right.end())
            {
                log_fatal("The PROPOSALParticle '%s' can not be converted to a I3Particle", particle_def.name.c_str());
            } else
            {
                return proposal_iterator->second;
            }
            break;
        }
        case DynamicData::Brems:
            return I3Particle::Brems;
            break;
        case DynamicData::Epair:
            return I3Particle::PairProd;
            break;
        case DynamicData::DeltaE:
            return I3Particle::DeltaE;
            break;
        case DynamicData::NuclInt:
            return I3Particle::NuclInt;
            break;
        default:
            log_fatal("PROPOSAL Particle can not be converted to a I3Particle");
    }
}

// ------------------------------------------------------------------------- //
PROPOSAL::ParticleDef I3PROPOSALParticleConverter::GeneratePROPOSALType(const I3Particle::ParticleType& ptype_I3)
{
    using namespace PROPOSAL;

    switch (ptype_I3)
    {
        case I3Particle::MuMinus:
            return MuMinusDef::Get();
            break;
        case I3Particle::MuPlus:
            return MuPlusDef::Get();
            break;
        case I3Particle::TauMinus:
            return TauMinusDef::Get();
            break;
        case I3Particle::TauPlus:
            return TauPlusDef::Get();
            break;
        case I3Particle::EMinus:
            return EMinusDef::Get();
            break;
        case I3Particle::EPlus:
            return EPlusDef::Get();
            break;
        default:
        {

            I3Particle i3particle;
            i3particle.SetType(ptype_I3);

            log_fatal("The I3Particle '%s' with type '%i' can not be converted to a PROPOSALParticle",
                      i3particle.GetTypeString().c_str(),
                      ptype_I3);
        }
    }
}
