
#include <boost/assign.hpp>

#include "PROPOSAL-icetray/Converter.h"

I3PROPOSALParticleConverter::I3PROPOSALParticleConverter()
{
    i3_to_proposal_.emplace(I3Particle::MuMinus, PROPOSAL::MuMinusDef::Get());
    i3_to_proposal_.emplace(I3Particle::MuPlus, PROPOSAL::MuPlusDef::Get());
    i3_to_proposal_.emplace(I3Particle::TauMinus, PROPOSAL::TauMinusDef::Get());
    i3_to_proposal_.emplace(I3Particle::TauPlus, PROPOSAL::TauPlusDef::Get());
    i3_to_proposal_.emplace(I3Particle::EMinus, PROPOSAL::EMinusDef::Get());
    i3_to_proposal_.emplace(I3Particle::EPlus, PROPOSAL::EPlusDef::Get());
    i3_to_proposal_.emplace(I3Particle::NuMuBar, PROPOSAL::NuMuBarDef::Get());
    i3_to_proposal_.emplace(I3Particle::NuMu, PROPOSAL::NuMuDef::Get());
    i3_to_proposal_.emplace(I3Particle::NuTauBar, PROPOSAL::NuTauBarDef::Get());
    i3_to_proposal_.emplace(I3Particle::NuTau, PROPOSAL::NuTauDef::Get());
    i3_to_proposal_.emplace(I3Particle::NuEBar, PROPOSAL::NuEBarDef::Get());
    i3_to_proposal_.emplace(I3Particle::NuE, PROPOSAL::NuEDef::Get());
    i3_to_proposal_.emplace(I3Particle::Monopole, PROPOSAL::MonopoleDef::Get());
    i3_to_proposal_.emplace(I3Particle::STauMinus, PROPOSAL::StauMinusDef::Get());
    i3_to_proposal_.emplace(I3Particle::STauPlus, PROPOSAL::StauPlusDef::Get());
    i3_to_proposal_.emplace(I3Particle::SMPMinus, PROPOSAL::SMPMinusDef::Get());
    i3_to_proposal_.emplace(I3Particle::SMPPlus, PROPOSAL::SMPPlusDef::Get());
    i3_to_proposal_.emplace(I3Particle::Pi0, PROPOSAL::Pi0Def::Get());
    i3_to_proposal_.emplace(I3Particle::PiMinus, PROPOSAL::PiMinusDef::Get());
    i3_to_proposal_.emplace(I3Particle::PiPlus, PROPOSAL::PiPlusDef::Get());
    i3_to_proposal_.emplace(I3Particle::K0_Short, PROPOSAL::K0Def::Get());
    i3_to_proposal_.emplace(I3Particle::K0_Long, PROPOSAL::K0Def::Get());
    i3_to_proposal_.emplace(I3Particle::KMinus, PROPOSAL::KMinusDef::Get());
    i3_to_proposal_.emplace(I3Particle::KPlus, PROPOSAL::KPlusDef::Get());
}

// ------------------------------------------------------------------------- //
I3Particle::ParticleType I3PROPOSALParticleConverter::GenerateI3Type(const PROPOSAL::DynamicData& secondary) const
{
    PROPOSAL::DynamicData::Type type = secondary.GetTypeId();

    switch (type)
    {
        case PROPOSAL::DynamicData::Particle:
        {
            const PROPOSAL::Particle& particle = static_cast<const PROPOSAL::Particle&>(secondary);
            PROPOSAL::ParticleDef particle_def = particle.GetParticleDef();

            for (const auto &pair : i3_to_proposal_) {
                if (particle_def == pair.second) {
                    return pair.first;
                }
            }
            log_warn("The PROPOSAL Particle '%s' can not be converted to a I3Particle", particle_def.name.c_str());
            return I3Particle::unknown;
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
const PROPOSAL::ParticleDef& I3PROPOSALParticleConverter::GeneratePROPOSALType(const I3Particle::ParticleType& ptype_I3) const
{
    auto it = i3_to_proposal_.find(ptype_I3);
    if (it == i3_to_proposal_.cend()) {
        I3Particle i3particle;
        i3particle.SetType(ptype_I3);

        log_warn("The I3Particle '%s' with type '%i' can not be converted to a PROPOSAL Particle",
                  i3particle.GetTypeString().c_str(),
                  ptype_I3);

        // Return an empty particle definition
        return null_definition_;
    } else {
        return it->second;
    }
}

// ------------------------------------------------------------------------- //
PROPOSAL::Particle I3PROPOSALParticleConverter::GeneratePROPOSALParticle(const I3Particle& p) const
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
    if (std::isnan(length))
    {
        length = 0.0;
    }
    else if (std::isinf(length))
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
I3Particle I3PROPOSALParticleConverter::GenerateI3Particle(const PROPOSAL::DynamicData& pp) const
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
