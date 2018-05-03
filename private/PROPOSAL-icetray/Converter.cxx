
#include <boost/assign.hpp>

#include "PROPOSAL-icetray/Converter.h"

namespace I3PROPOSALParticleConverter {

typedef boost::bimap<I3Particle::ParticleType, std::string> bimap_ParticleType;
static const bimap_ParticleType I3_PROPOSAL_ParticleType_bimap =
    boost::assign::list_of<bimap_ParticleType::relation>
        (I3Particle::MuMinus, "MuMinus")
        (I3Particle::MuPlus, "MuPlus")
        (I3Particle::TauMinus, "TauMinus")
        (I3Particle::TauPlus, "TauPlus")
        (I3Particle::EMinus, "EMinus")
        (I3Particle::EPlus, "EPlus")
        (I3Particle::NuMu, "NuMu")
        (I3Particle::NuMuBar, "NuMuBar")
        (I3Particle::NuE, "NuE")
        (I3Particle::NuEBar, "NuEBar")
        (I3Particle::NuTau, "NuTau")
        (I3Particle::NuTauBar, "NuTauBar")
        (I3Particle::Brems, "Brems")
        (I3Particle::DeltaE, "DeltaE")
        (I3Particle::PairProd, "EPair")
        (I3Particle::NuclInt, "NuclInt")
        (I3Particle::MuPair, "MuPair")
        (I3Particle::Hadrons, "Hadrons")
        (I3Particle::Monopole, "Monopole")
        (I3Particle::STauMinus, "STauMinus")
        (I3Particle::STauPlus, "STauPlus")
        (I3Particle::Gamma, "Gamma")
        (I3Particle::Pi0, "Pi0")
        (I3Particle::PiPlus, "PiPlus")
        (I3Particle::PiMinus, "PiMinus")
        (I3Particle::K0_Short, "K0") // TODO(mario):  Fri 2017/11/17
        (I3Particle::KPlus, "KPlus")
        (I3Particle::KMinus, "KMinus")
        (I3Particle::PPlus, "PPlus")
        (I3Particle::PMinus, "PMinus");

} // namespace I3PROPOSALParticleConverter

// ------------------------------------------------------------------------- //
I3Particle::ParticleType I3PROPOSALParticleConverter::GenerateI3Type(const PROPOSAL::DynamicData& secondary)
{
    PROPOSAL::DynamicData::Type type = secondary.GetTypeId();

    switch (type)
    {
        case PROPOSAL::DynamicData::Particle:
        {
            const PROPOSAL::Particle& particle = static_cast<const PROPOSAL::Particle&>(secondary);
            PROPOSAL::ParticleDef particle_def = particle.GetParticleDef();

            // I3Particle::ParticleType ptype_I3;

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
        case PROPOSAL::DynamicData::Brems:
            return I3Particle::Brems;
            break;
        case PROPOSAL::DynamicData::Epair:
            return I3Particle::PairProd;
            break;
        case PROPOSAL::DynamicData::DeltaE:
            return I3Particle::DeltaE;
            break;
        case PROPOSAL::DynamicData::NuclInt:
            return I3Particle::NuclInt;
            break;
        case PROPOSAL::DynamicData::ContinuousEnergyLoss:
            return I3Particle::ContinuousEnergyLoss;
            break;
        default:
            log_fatal("PROPOSAL Particle can not be converted to a I3Particle");
    }
}

// ------------------------------------------------------------------------- //
PROPOSAL::ParticleDef I3PROPOSALParticleConverter::GeneratePROPOSALType(const I3Particle::ParticleType& ptype_I3)
{
    switch (ptype_I3)
    {
        case I3Particle::MuMinus:
            return PROPOSAL::MuMinusDef::Get();
            break;
        case I3Particle::MuPlus:
            return PROPOSAL::MuPlusDef::Get();
            break;
        case I3Particle::TauMinus:
            return PROPOSAL::TauMinusDef::Get();
            break;
        case I3Particle::TauPlus:
            return PROPOSAL::TauPlusDef::Get();
            break;
        case I3Particle::EMinus:
            return PROPOSAL::EMinusDef::Get();
            break;
        case I3Particle::EPlus:
            return PROPOSAL::EPlusDef::Get();
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

// ------------------------------------------------------------------------- //
PROPOSAL::Particle I3PROPOSALParticleConverter::GeneratePROPOSALParticle(const I3Particle& p)
{
    /**
     * Natural units of PROPOSAL is cm, rad, MeV, and s.
     * Therefore we need to convert explicitly to
     * PROPOSAL units before passing the propagate method
     */

    double x = p.GetPos().GetX() / I3Units::cm;
    double y = p.GetPos().GetY() / I3Units::cm;
    double z = p.GetPos().GetZ() / I3Units::cm;
    double theta = p.GetDir().CalcTheta() / I3Units::radian;
    double phi = p.GetDir().CalcPhi() / I3Units::radian;

    // log_debug("Name of particle to propagate: %s", PROPOSAL::Particle::GetName(GeneratePROPOSALType(p)).c_str());

    PROPOSAL::ParticleDef particle_def = I3PROPOSALParticleConverter::GeneratePROPOSALType(p.GetType());
    PROPOSAL::Particle particle(particle_def);
    particle.SetPosition(PROPOSAL::Vector3D(x, y, z));

    PROPOSAL::Vector3D direction;
    direction.SetSphericalCoordinates(1.0, phi, theta);
    direction.CalculateCartesianFromSpherical();
    particle.SetDirection(direction);

    particle.SetEnergy(p.GetEnergy() / I3Units::MeV);
    particle.SetTime(p.GetTime() / I3Units::second);
    particle.SetPropagatedDistance(p.GetLength() / I3Units::cm);

    return particle;
}

// ------------------------------------------------------------------------- //
I3Particle I3PROPOSALParticleConverter::GenerateI3Particle(const PROPOSAL::DynamicData& pp)
{
    double x = pp.GetPosition().GetX() * I3Units::cm;
    double y = pp.GetPosition().GetY() * I3Units::cm;
    double z = pp.GetPosition().GetZ() * I3Units::cm;

    double theta = pp.GetDirection().GetTheta() * I3Units::radian;
    double phi   = pp.GetDirection().GetPhi() * I3Units::radian;

    I3Particle i3_particle;
    i3_particle.SetType(I3PROPOSALParticleConverter::GenerateI3Type(pp));
    i3_particle.SetLocationType(I3Particle::InIce);

    i3_particle.SetPos(x, y, z);
    i3_particle.SetThetaPhi(theta, phi);

    i3_particle.SetLength(pp.GetPropagatedDistance() * I3Units::cm);
    i3_particle.SetTime(pp.GetTime() * I3Units::second);
    i3_particle.SetEnergy(pp.GetEnergy() * I3Units::MeV);
    
    log_trace("MMC DEBUG SEC \n  pos=(%g,%g,%g) ang=(%g,%g)  e=%g t=%g  l=%g",
        x, y, z, theta, phi,
        pp.GetEnergy() * I3Units::MeV,
        pp.GetTime() * I3Units::second,
        pp.GetPropagatedDistance() * I3Units::cm);

    return i3_particle;
}
