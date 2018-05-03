/** $Id: MuonPropagator.cxx 137064 2015-08-31 18:24:47Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 137064 $
 * $Date: 2015-08-31 20:24:47 +0200 (Mo, 31. Aug 2015) $
 */

#include <boost/assign.hpp>
#include <boost/foreach.hpp>

#include "MuonGun/MuonPropagator.h"

namespace I3MuonGun {

MuonPropagator::MuonPropagator(const std::string& medium, double ecut, double vcut, double rho)
    : propagator_(NULL)
{
    PROPOSAL::Sector::Definition sec_def;

    sec_def.cut_settings.SetEcut(ecut);
    sec_def.cut_settings.SetVcut(vcut);

    sec_def.SetMedium(*PROPOSAL::MediumFactory::Get().CreateMedium(medium, 1.0));

    PROPOSAL::Sphere detector(PROPOSAL::Vector3D(0.0, 0.0, 0.0), 0.0, 1e18);
    sec_def.SetGeometry(detector);

    std::vector<PROPOSAL::Sector::Definition> sector_definitions;
    sector_definitions.push_back(sec_def);

    PROPOSAL::InterpolationDef interpol_def;

    std::ostringstream prefix;
    prefix << getenv("I3_BUILD") << "/MuonGun/resources/tables/icecube";

    interpol_def.path_to_tables = prefix.str();

    propagator_ = new PROPOSAL::Propagator(PROPOSAL::MuMinusDef::Get(),
                                           sector_definitions,
                                           detector,
                                           interpol_def);
}

MuonPropagator::~MuonPropagator()
{
    delete propagator_;
}

void MuonPropagator::SetSeed(int seed)
{
    PROPOSAL::RandomGenerator::Get().SetSeed(seed);
}

inline std::string GetMMCName(I3Particle::ParticleType pt)
{
    std::string name;

    if (pt == I3Particle::MuMinus)
        name = "mu-";
    else if (pt == I3Particle::MuPlus)
        name = "mu+";

    return name;
}

std::string MuonPropagator::GetName(const I3Particle& p)
{
    return GetMMCName(p.GetType());
}


/** Differential stochastic rate: d^2N/dv/dx [1/m] */
double MuonPropagator::GetStochasticRate(double energy, double fraction, I3Particle::ParticleType type) const
{
    // the propagator_ is always initialised with a Muon
    // there is no need to specify the type of the particle
    // this input parameter is historic and can be vanished
    // it just exist for backward compatibility
    (void)type;

    // Check kinematics
    if (fraction <= 0 || energy*(1-fraction) <= propagator_->GetParticle().GetMass()*I3Units::MeV)
        return 0.;

    double contrib;
    double rate = 0.0;

    const std::vector<PROPOSAL::CrossSection*>& crosssections = propagator_->GetCurrentSector()->GetUtility().GetCrosssections();

    for (std::vector<PROPOSAL::CrossSection*>::const_iterator iter = crosssections.begin(); iter != crosssections.end(); ++iter)
    {
        if((*iter)->GetTypeId() == PROPOSAL::DynamicData::Brems || 
            (*iter)->GetTypeId() == PROPOSAL::DynamicData::Epair || 
            (*iter)->GetTypeId() == PROPOSAL::DynamicData::NuclInt)
        {
            for (int i = 0; i < propagator_->GetCurrentSector()->GetUtility().GetMedium().GetNumComponents(); ++i)
            {
                (*iter)->GetParametrization().SetCurrentComponent(i);
                contrib = (*iter)->GetParametrization().DifferentialCrossSection(energy, fraction);
                if(std::isfinite(contrib) && contrib > 0)
                    rate += contrib;
            }
        }
        else // if((*iter)->GetTypeId() == PROPOSAL::DynamicData::DeltaE)
        {
            contrib = (*iter)->GetParametrization().DifferentialCrossSection(energy, fraction);
            if(std::isfinite(contrib) && contrib > 0)
                rate += contrib;
        }
    }
    return rate / I3Units::cm;
}

/** total stochastic rate: dN/dx [1/m] */
double MuonPropagator::GetTotalStochasticRate(double energy, I3Particle::ParticleType type) const
{
    // the propagator_ is always initialised with a Muon
    // there is no need to specify the type of the particle
    // this input parameter is historic and can be vanished
    // it just exist for backward compatibility
    (void)type;

    double rate = 0.0;
    const std::vector<PROPOSAL::CrossSection*>& crosssections = propagator_->GetCurrentSector()->GetUtility().GetCrosssections();

    for (std::vector<PROPOSAL::CrossSection*>::const_iterator iter = crosssections.begin(); iter != crosssections.end(); ++iter)
    {
        rate += (*iter)->CalculatedNdx(energy);
    }
    return rate / I3Units::cm;
}

I3Particle MuonPropagator::propagate(const I3Particle& p,
                                     double distance,
                                     boost::shared_ptr<std::vector<I3Particle> > losses)
{
    I3Particle endpoint(p);

    double x, y, z, theta, phi;

    PROPOSAL::Particle pp = propagator_->GetParticle();
    pp.SetParentParticleId(0);
    pp.SetParticleId(0);

    x = p.GetPos().GetX() / I3Units::cm;
    y = p.GetPos().GetY() / I3Units::cm;
    z = p.GetPos().GetZ() / I3Units::cm;

    theta = p.GetDir().CalcTheta() / I3Units::radian;
    phi   = p.GetDir().CalcPhi() / I3Units::radian;

    pp.SetPosition(PROPOSAL::Vector3D(x, y, z));

    PROPOSAL::Vector3D direction;
    direction.SetSphericalCoordinates(1.0, phi, theta);
    direction.CalculateCartesianFromSpherical();
    pp.SetDirection(direction);

    pp.SetEnergy(p.GetEnergy() / I3Units::MeV);
    pp.SetTime(p.GetTime() / I3Units::second);
    pp.SetPropagatedDistance(p.GetLength() / I3Units::cm);

    propagator_->Propagate(distance / I3Units::cm);

    x = pp.GetPosition().GetX() * I3Units::cm;
    y = pp.GetPosition().GetY() * I3Units::cm;
    z = pp.GetPosition().GetZ() * I3Units::cm;

    theta = pp.GetDirection().GetTheta() * I3Units::radian;
    phi   = pp.GetDirection().GetPhi() * I3Units::radian;

    endpoint.SetPos(x, y, z);
    endpoint.SetThetaPhi(theta, phi);

    endpoint.SetEnergy(pp.GetEnergy() * I3Units::MeV);
    endpoint.SetLength(pp.GetPropagatedDistance() * I3Units::cm);
    endpoint.SetTime(pp.GetTime() * I3Units::second);

    if (losses)
    {
        std::vector<PROPOSAL::DynamicData*> history = PROPOSAL::Output::getInstance().GetSecondarys();
        for (unsigned int i = 0; i < history.size(); i++)
        {
            losses->push_back(GenerateI3Particle(*history[i]));
        }
    }

    PROPOSAL::Output::getInstance().ClearSecondaryVector();

    return endpoint;
}

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

// ------------------------------------------------------------------------- //
I3Particle::ParticleType MuonPropagator::GenerateI3Type(const PROPOSAL::DynamicData& secondary)
{
    PROPOSAL::DynamicData::Type type = secondary.GetTypeId();

    switch (type)
    {
        case PROPOSAL::DynamicData::Particle:
        {
            const PROPOSAL::Particle& particle = static_cast<const PROPOSAL::Particle&>(secondary);
            PROPOSAL::ParticleDef particle_def = particle.GetParticleDef();

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
I3Particle MuonPropagator::GenerateI3Particle(const PROPOSAL::DynamicData& pp)
{
    double x = pp.GetPosition().GetX() * I3Units::cm;
    double y = pp.GetPosition().GetY() * I3Units::cm;
    double z = pp.GetPosition().GetZ() * I3Units::cm;

    double theta = pp.GetDirection().GetTheta() * I3Units::radian;
    double phi   = pp.GetDirection().GetPhi() * I3Units::radian;

    I3Particle i3_particle;
    i3_particle.SetType(GenerateI3Type(pp));
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


void Crust::AddLayer(I3Surfaces::SurfacePtr s, boost::shared_ptr<MuonPropagator> p)
{
    boundaries_.push_back(s);
    propagators_.push_back(p);
}

I3Particle Crust::Ingest(const I3Particle& p)
{
    I3Particle propped(p);
    double l = 0;
    for (unsigned i = 0; (propped.GetEnergy() > 0) && (i < boundaries_.size()); i++)
    {
        double dx = boundaries_[i]->GetIntersection(propped.GetPos(), propped.GetDir()).first;
        if (dx > 0)
            propped = (i > 0 ? propagators_[i - 1] : defaultPropagator_)->propagate(propped, dx);
        // Force lengths to measure the distance back to the outermost surface
        if (i > 0)
            l += std::min(dx, propped.GetLength());
    }
    propped.SetLength(l);

    return propped;
}

} // namespace I3MuonGun
