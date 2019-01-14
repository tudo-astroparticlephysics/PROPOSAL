
#include <boost/assign.hpp>

#include "PROPOSAL-icetray/Converter.h"


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

            if (particle_def == PROPOSAL::MuMinusDef::Get())
                return I3Particle::MuMinus;
            else if (particle_def == PROPOSAL::MuPlusDef::Get())
                return I3Particle::MuPlus;
            else if (particle_def == PROPOSAL::MuPlusDef::Get())
                return I3Particle::MuPlus;
            else if (particle_def == PROPOSAL::TauPlusDef::Get())
                return I3Particle::TauPlus;
            else if (particle_def == PROPOSAL::TauMinusDef::Get())
                return I3Particle::TauMinus;
            else if (particle_def == PROPOSAL::EPlusDef::Get())
                return I3Particle::EPlus;
            else if (particle_def == PROPOSAL::EMinusDef::Get())
                return I3Particle::EMinus;
            else if (particle_def == PROPOSAL::NuMuBarDef::Get())
                return I3Particle::NuMuBar;
            else if (particle_def == PROPOSAL::NuMuDef::Get())
                return I3Particle::NuMu;
            else if (particle_def == PROPOSAL::NuTauBarDef::Get())
                return I3Particle::NuTauBar;
            else if (particle_def == PROPOSAL::NuTauDef::Get())
                return I3Particle::NuTau;
            else if (particle_def == PROPOSAL::NuEBarDef::Get())
                return I3Particle::NuEBar;
            else if (particle_def == PROPOSAL::NuEDef::Get())
                return I3Particle::NuE;
            else if (particle_def == PROPOSAL::MonopoleDef::Get())
                return I3Particle::Monopole;
            else if (particle_def == PROPOSAL::StauPlusDef::Get())
                return I3Particle::STauPlus;
            else if (particle_def == PROPOSAL::StauMinusDef::Get())
                return I3Particle::STauMinus;
            else if (particle_def == PROPOSAL::SMPPlusDef::Get())
                return I3Particle::SMPPlus;
            else if (particle_def == PROPOSAL::SMPMinusDef::Get())
                return I3Particle::SMPMinus;
            else if (particle_def == PROPOSAL::Pi0Def::Get())
                return I3Particle::Pi0;
            else if (particle_def == PROPOSAL::PiPlusDef::Get())
                return I3Particle::PiPlus;
            else if (particle_def == PROPOSAL::PiMinusDef::Get())
                return I3Particle::PiMinus;
            else if (particle_def == PROPOSAL::K0Def::Get())
                return I3Particle::K0_Short;
            else if (particle_def == PROPOSAL::KPlusDef::Get())
                return I3Particle::KPlus;
            else if (particle_def == PROPOSAL::KMinusDef::Get())
                return I3Particle::KMinus;
            else
            {
                log_warn("The PROPOSAL Particle '%s' can not be converted to a I3Particle", particle_def.name.c_str());
                return I3Particle::unknown;
            }
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
        case PROPOSAL::DynamicData::MuPair:
            return I3Particle::MuPair;
            break;
        case PROPOSAL::DynamicData::ContinuousEnergyLoss:
            return I3Particle::ContinuousEnergyLoss;
            break;
        default:
            {
                log_warn("PROPOSAL Particle can not be converted to a I3Particle");
                return I3Particle::unknown;
            }
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
        case I3Particle::NuMuBar:
            return PROPOSAL::NuMuBarDef::Get();
            break;
        case I3Particle::NuMu:
            return PROPOSAL::NuMuDef::Get();
            break;
        case I3Particle::NuTauBar:
            return PROPOSAL::NuTauBarDef::Get();
            break;
        case I3Particle::NuTau:
            return PROPOSAL::NuTauDef::Get();
            break;
        case I3Particle::NuEBar:
            return PROPOSAL::NuEBarDef::Get();
            break;
        case I3Particle::NuE:
            return PROPOSAL::NuEDef::Get();
            break;
        case I3Particle::Monopole:
            return PROPOSAL::MonopoleDef::Get();
            break;
        case I3Particle::STauMinus:
            return PROPOSAL::StauMinusDef::Get();
            break;
        case I3Particle::STauPlus:
            return PROPOSAL::StauPlusDef::Get();
            break;
        case I3Particle::SMPMinus:
            return PROPOSAL::SMPMinusDef::Get();
            break;
        case I3Particle::SMPPlus:
            return PROPOSAL::SMPPlusDef::Get();
            break;
        case I3Particle::Pi0:
            return PROPOSAL::Pi0Def::Get();
            break;
        case I3Particle::PiMinus:
            return PROPOSAL::PiMinusDef::Get();
            break;
        case I3Particle::PiPlus:
            return PROPOSAL::PiPlusDef::Get();
            break;
        case I3Particle::K0_Long:
        case I3Particle::K0_Short:
            return PROPOSAL::K0Def::Get();
            break;
        case I3Particle::KMinus:
            return PROPOSAL::KMinusDef::Get();
            break;
        case I3Particle::KPlus:
            return PROPOSAL::KPlusDef::Get();
            break;
        default:
        {
            I3Particle i3particle;
            i3particle.SetType(ptype_I3);

            log_warn("The I3Particle '%s' with type '%i' can not be converted to a PROPOSAL Particle",
                      i3particle.GetTypeString().c_str(),
                      ptype_I3);

            // Return an empty particle definition
            return PROPOSAL::ParticleDef::Builder().build();
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
    double phi   = p.GetDir().CalcPhi() / I3Units::radian;

    double energy = p.GetEnergy() / I3Units::MeV;
    double time   = p.GetTime() / I3Units::second;
    double length = p.GetLength() / I3Units::cm;

    // The Muons from NuGen have NaN as default propagated length.
    // So this has to be corrected.
    if (isnan(length))
    {
        length = 0.0;
    }
    else if (isinf(length))
    {
        log_fatal("the propagated length is Inf, should be finite or NaN.");
    }

    // log_debug("Name of particle to propagate: %s", PROPOSAL::Particle::GetName(GeneratePROPOSALType(p)).c_str());

    PROPOSAL::ParticleDef particle_def = I3PROPOSALParticleConverter::GeneratePROPOSALType(p.GetType());
    PROPOSAL::Particle particle(particle_def);
    particle.SetPosition(PROPOSAL::Vector3D(x, y, z));

    PROPOSAL::Vector3D direction;
    direction.SetSphericalCoordinates(1.0, phi, theta);
    direction.CalculateCartesianFromSpherical();
    particle.SetDirection(direction);

    particle.SetEnergy(energy);
    particle.SetTime(time);
    particle.SetPropagatedDistance(length);

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

    double energy = pp.GetEnergy() * I3Units::MeV;
    double time   = pp.GetTime() * I3Units::second;
    double length = pp.GetPropagatedDistance() * I3Units::cm;

    I3Particle i3_particle;
    i3_particle.SetType(I3PROPOSALParticleConverter::GenerateI3Type(pp));
    i3_particle.SetLocationType(I3Particle::InIce);

    i3_particle.SetPos(x, y, z);
    i3_particle.SetThetaPhi(theta, phi);

    i3_particle.SetEnergy(energy);
    i3_particle.SetTime(time);
    i3_particle.SetLength(length);

    log_trace("MMC DEBUG SEC \n  pos=(%g,%g,%g) ang=(%g,%g)  e=%g t=%g  l=%g",
        x, y, z, theta, phi,
        energy, time, length);

    return i3_particle;
}
