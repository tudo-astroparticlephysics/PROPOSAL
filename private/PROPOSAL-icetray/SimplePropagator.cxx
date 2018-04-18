/** $Id: SimplePropagator.cxx 150208 2016-09-21 11:32:47Z tfuchs $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 150208 $
 * $Date: 2016-09-21 13:32:47 +0200 (Mi, 21. Sep 2016) $
 */

#include <boost/assign.hpp>
#include <boost/bimap.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>

#include "PROPOSAL-icetray/Converter.h"
#include "PROPOSAL-icetray/SimplePropagator.h"

#include "PROPOSAL/Output.h"

using namespace PROPOSAL;

// ParticleType::Enum GetPROPOSALType(I3Particle::ParticleType pt)
// {
//     ParticleType::Enum particle_type;
//     switch (pt) {
//         case I3Particle::MuMinus:
//             particle_type = ParticleType::MuMinus;
//             break;
//         case I3Particle::MuPlus:
//             particle_type = ParticleType::MuPlus;
//             break;
//         case I3Particle::TauMinus:
//             particle_type = ParticleType::TauMinus;
//             break;
//         case I3Particle::TauPlus:
//             particle_type = ParticleType::TauPlus;
//             break;
//         default:
//             log_fatal_stream("Unsupported particle type: " << pt);
//     }
//     return particle_type;
// }

SimplePropagator::SimplePropagator(I3Particle::ParticleType pt,
                                   const std::string& medium,
                                   double ecut,
                                   double vcut,
                                   double rho)
{
    std::ostringstream prefix;
    prefix << getenv("I3_BUILD") << "/MuonGun/resources/tables/icecube";

    // Sector definition

    PROPOSAL::SectorFactory::Definition sec_def;

    sec_def.stopping_decay              = true;
    sec_def.scattering_model            = PROPOSAL::ScatteringFactory::Moliere;
    sec_def.do_exact_time_calculation   = true;
    sec_def.do_continuous_randomization = ecut < 0;

    // Medium

    sec_def.medium_def.type               = PROPOSAL::MediumFactory::Get().GetEnumFromString(medium);
    sec_def.medium_def.density_correction = rho;

    // Geometry

    PROPOSAL::Vector3D position(0.0, 0.0, 0.0);

    PROPOSAL::Sphere detector(position, 0.0, 1e18);

    sec_def.geometry_def.shape        = PROPOSAL::GeometryFactory::Sphere;
    sec_def.geometry_def.position     = position;
    sec_def.geometry_def.inner_radius = 0.0;
    sec_def.geometry_def.radius       = 1e18;

    // Cuts

    sec_def.e_cut = ecut;
    sec_def.v_cut = vcut;

    // Parametrizations

    sec_def.utility_def.brems_def.parametrization = PROPOSAL::BremsstrahlungFactory::KelnerKokoulinPetrukhin;

    sec_def.utility_def.photo_def.parametrization = PROPOSAL::PhotonuclearFactory::AbramowiczLevinLevyMaor97;

    sec_def.utility_def.epair_def.lpm_effect = true;
    sec_def.utility_def.brems_def.lpm_effect = true;

    // Interpolation

    PROPOSAL::InterpolationDef interpolation_def;
    interpolation_def.path_to_tables = prefix.str();

    // Init new propagator

    // TODO(mario): Check for muon, tau Mon 2017/11/06
    propagator_ = new PROPOSAL::Propagator(PROPOSAL::MuMinusDef::Get(),
                                           boost::assign::list_of<PROPOSAL::SectorFactory::Definition>(sec_def),
                                           detector,
                                           interpolation_def);
}

SimplePropagator::~SimplePropagator()
{
    delete propagator_;
}

void SimplePropagator::SetSeed(int seed)
{
    PROPOSAL::RandomGenerator::Get().SetSeed(seed);
}

void SimplePropagator::SetRandomNumberGenerator(I3RandomServicePtr rng)
{
    boost::function<double()> f = boost::bind(&I3RandomService::Uniform, rng, 0, 1);
    PROPOSAL::RandomGenerator::Get().SetRandomNumberGenerator(f);
}

// typedef boost::bimap<I3Particle::ParticleType, std::string> bimap_ParticleType;
// static const bimap_ParticleType I3_PROPOSAL_ParticleType_bimap = boost::assign::list_of<bimap_ParticleType::relation>
//     (I3Particle::MuMinus,   "MuMinus")
//     (I3Particle::MuPlus,    "MuPlus")
//     (I3Particle::TauMinus,  "TauMinus")
//     (I3Particle::TauPlus,   "TauPlus")
//     (I3Particle::EMinus,    "EMinus")
//     (I3Particle::EPlus,     "EPlus")
//     (I3Particle::NuMu,      "NuMu")
//     (I3Particle::NuMuBar,   "NuMuBar")
//     (I3Particle::NuE,       "NuE")
//     (I3Particle::NuEBar,    "NuEBar")
//     (I3Particle::NuTau,     "NuTau")
//     (I3Particle::NuTauBar,  "NuTauBar")
//     (I3Particle::Brems,     "Brems")
//     (I3Particle::DeltaE,    "DeltaE")
//     (I3Particle::PairProd,  "EPair")
//     (I3Particle::NuclInt,   "NuclInt")
//     (I3Particle::MuPair,    "MuPair")
//     (I3Particle::Hadrons,   "Hadrons")
//     (I3Particle::Monopole,  "Monopole")
//     (I3Particle::STauMinus, "STauMinus")
//     (I3Particle::STauPlus,  "STauPlus")
//     (I3Particle::Gamma,     "Gamma")
//     (I3Particle::Pi0,       "Pi0")
//     (I3Particle::PiPlus,    "PiPlus")
//     (I3Particle::PiMinus,   "PiMinus")
//     (I3Particle::KPlus,     "KPlus")
//     (I3Particle::KMinus,    "KMinus")
//     (I3Particle::PPlus,     "PPlus")
//     (I3Particle::PMinus,    "PMinus");

// ------------------------------------------------------------------------- //
// ParticleType::Enum SimplePropagator::GeneratePROPOSALType(const I3Particle& p)
// {
//     I3Particle::ParticleType ptype_I3 = p.GetType();
//     ParticleType::Enum ptype_PROPOSAL;
//
//     bimap_ParticleType::left_const_iterator i3_iterator = I3_PROPOSAL_ParticleType_bimap.left.find(ptype_I3);
//     if (i3_iterator == I3_PROPOSAL_ParticleType_bimap.left.end())
//     {
//         log_fatal("The I3Particle '%s' with type '%i' can not be converted to a PROPOSALParticle"
//             , p.GetTypeString().c_str(), ptype_I3);
//     }
//
//     ptype_PROPOSAL = I3_PROPOSAL_ParticleType_bimap.left.find(ptype_I3) -> second;
//
//     return ptype_PROPOSAL;
// }

// I3Particle::ParticleType SimplePropagator::GenerateI3Type(const PROPOSAL::DynamicData& secondary)
// {
//     DynamicData::Type type = secondary.GetTypeId();
//
//     if (type == DynamicData::Particle)
//     {
//         const PROPOSAL::Particle& particle = static_cast<const PROPOSAL::Particle&>(secondary);
//         ParticleDef particle_def = particle.GetParticleDef();
//
//         I3Particle::ParticleType ptype_I3;
//
//         bimap_ParticleType::right_const_iterator proposal_iterator =
//         I3_PROPOSAL_ParticleType_bimap.right.find(particle_def.name); if (proposal_iterator ==
//         I3_PROPOSAL_ParticleType_bimap.right.end())
//         {
//             log_fatal("The PROPOSALParticle '%s' can not be converted to a I3Particle", particle_def.name.c_str());
//         }
//         else
//         {
//             return proposal_iterator->second;
//         }
//     }
//     else if(type == DynamicData::Brems) return I3Particle::Brems;
//     else if(type == DynamicData::Epair) return I3Particle::PairProd;
//     else if(type == DynamicData::DeltaE) return I3Particle::DeltaE;
//     else if(type == DynamicData::NuclInt) return I3Particle::NuclInt;
//     else
//     {
//         log_fatal("PROPOSAL Particle can not be converted to a I3Particle");
//     }
// }

I3Particle SimplePropagator::to_I3Particle(const PROPOSAL::DynamicData& pp)
{
    I3Particle p;

    double x = pp.GetPosition().GetX() * I3Units::cm;
    double y = pp.GetPosition().GetY() * I3Units::cm;
    double z = pp.GetPosition().GetZ() * I3Units::cm;

    double theta = pp.GetDirection().GetTheta() * I3Units::degree;
    double phi   = pp.GetDirection().GetPhi() * I3Units::degree;

    p.SetPos(x, y, z);
    p.SetThetaPhi(theta, phi);
    p.SetLength(pp.GetPropagatedDistance() * I3Units::cm);
    p.SetTime(pp.GetTime() * I3Units::second);

    p.SetType(I3PROPOSALParticleConverter::GenerateI3Type(pp));
    p.SetLocationType(I3Particle::InIce);
    p.SetEnergy(pp.GetEnergy() * I3Units::MeV);

    return p;
}

I3Particle SimplePropagator::propagate(const I3Particle& p,
                                       double distance,
                                       boost::shared_ptr<std::vector<I3Particle> > losses)
{
    I3Particle endpoint(p);

    double x, y, z, theta, phi = 0.0;

    PROPOSAL::Particle pp = propagator_->GetParticle();
    pp.SetParentParticleId(0);
    pp.SetParticleId(0);
    pp.SetTime(p.GetTime() / I3Units::second);

    x = p.GetPos().GetX() / I3Units::cm;
    y = p.GetPos().GetY() / I3Units::cm;
    z = p.GetPos().GetZ() / I3Units::cm;

    pp.SetPosition(PROPOSAL::Vector3D(x, y, z));
    pp.SetEnergy(p.GetEnergy() / I3Units::MeV);

    propagator_->Propagate(distance / I3Units::cm);

    endpoint.SetEnergy(pp.GetEnergy() * I3Units::MeV);

    x = pp.GetPosition().GetX() * I3Units::cm;
    y = pp.GetPosition().GetY() * I3Units::cm;
    z = pp.GetPosition().GetZ() * I3Units::cm;

    theta = pp.GetDirection().GetTheta() * I3Units::degree;
    phi   = pp.GetDirection().GetPhi() * I3Units::degree;

    endpoint.SetPos(x, y, z);
    endpoint.SetThetaPhi(theta, phi);
    endpoint.SetLength(pp.GetPropagatedDistance() * I3Units::cm);
    endpoint.SetTime(pp.GetTime() * I3Units::second);

    // Tomasz
    if (losses)
    {
        std::vector<PROPOSAL::DynamicData*> history = Output::getInstance().GetSecondarys();
        for (unsigned int i = 0; i < history.size(); i++)
        {
            losses->push_back(to_I3Particle(*history[i]));
        }
    }

    Output::getInstance().ClearSecondaryVector();

    return endpoint;
}
