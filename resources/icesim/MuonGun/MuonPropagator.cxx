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

MuonPropagator::MuonPropagator(const std::string &medium, double ecut, double vcut, double rho)
	: propagator_(NULL)
{
    PROPOSAL::SectorFactory::Definition sector_def;
    sector_def.e_cut = ecut;
    sector_def.v_cut = vcut;
    sector_def.medium_def.type = PROPOSAL::MediumFactory::Get().GetEnumFromString(medium);
    sector_def.medium_def.density_correction = rho;

    std::vector<PROPOSAL::SectorFactory::Definition> sector_definitions;
    sector_definitions.push_back(sector_def);

    PROPOSAL::InterpolationDef interpol_def;

	std::ostringstream prefix;
	prefix << getenv("I3_BUILD") << "/MuonGun/resources/tables/icecube";

    interpol_def.path_to_tables = prefix.str();

    propagator_ = new PROPOSAL::Propagator(PROPOSAL::MuMinusDef::Get(), sector_definitions, PROPOSAL::Sphere(PROPOSAL::Vector3D(), 1e18, 0.0), interpol_def);
}

MuonPropagator::~MuonPropagator()
{
	delete propagator_;
}

void
MuonPropagator::SetSeed(int seed)
{
    PROPOSAL::RandomGenerator::Get().SetSeed(seed);
}

inline std::string
GetMMCName(I3Particle::ParticleType pt)
{
	std::string name;

	if (pt == I3Particle::MuMinus)
		name="mu-";
	else if (pt == I3Particle::MuPlus)
		name="mu+";

	return name;
}

std::string
MuonPropagator::GetName(const I3Particle &p)
{
	return GetMMCName(p.GetType());
}

typedef boost::bimap<I3Particle::ParticleType, std::string> bimap_ParticleType;

static const bimap_ParticleType I3_PROPOSAL_ParticleType_bimap =
boost::assign::list_of<bimap_ParticleType::relation>
    (I3Particle::MuMinus,   "MuMinus")
    (I3Particle::MuPlus,    "MuPlus")
    (I3Particle::TauMinus,  "TauMinus")
    (I3Particle::TauPlus,   "TauPlus")
    (I3Particle::EMinus,    "EMinus")
    (I3Particle::EPlus,     "EPlus")
    (I3Particle::NuMu,      "NuMu")
    (I3Particle::NuMuBar,   "NuMuBar")
    (I3Particle::NuE,       "NuE")
    (I3Particle::NuEBar,    "NuEBar")
    (I3Particle::NuTau,     "NuTau")
    (I3Particle::NuTauBar,  "NuTauBar")
    (I3Particle::Brems,     "Brems")
    (I3Particle::DeltaE,    "DeltaE")
    (I3Particle::PairProd,  "EPair")
    (I3Particle::NuclInt,   "NuclInt")
    (I3Particle::MuPair,    "MuPair")
    (I3Particle::Hadrons,   "Hadrons")
    (I3Particle::Monopole,  "Monopole")
    (I3Particle::STauMinus, "STauMinus")
    (I3Particle::STauPlus,  "STauPlus")
    (I3Particle::Gamma,     "Gamma")
    (I3Particle::Pi0,       "Pi0")
    (I3Particle::PiPlus,    "PiPlus")
    (I3Particle::PiMinus,   "PiMinus")
    (I3Particle::KPlus,     "KPlus")
    (I3Particle::KMinus,    "KMinus")
    (I3Particle::PPlus,     "PPlus")
    (I3Particle::PMinus,    "PMinus");

// typedef std::map<int, I3Particle::ParticleType> particle_type_conversion_t;
// static const particle_type_conversion_t fromRDMCTable =
// boost::assign::list_of<std::pair<int, I3Particle::ParticleType> >
// (-100, I3Particle::unknown)
// (1, I3Particle::Gamma)
// (2, I3Particle::EPlus)
// (3, I3Particle::EMinus)
// (4, I3Particle::Nu)
// (5, I3Particle::MuPlus)
// (6, I3Particle::MuMinus)
// (7, I3Particle::Pi0)
// (8, I3Particle::PiPlus)
// (9, I3Particle::PiMinus)
// (11, I3Particle::KPlus)
// (12, I3Particle::KMinus)
// (14, I3Particle::PPlus)
// (15, I3Particle::PMinus)
// (33, I3Particle::TauPlus)
// (34, I3Particle::TauMinus)
// (41, I3Particle::Monopole)
// (201, I3Particle::NuE)
// (202, I3Particle::NuMu)
// (203, I3Particle::NuTau)
// (204, I3Particle::NuEBar)
// (205, I3Particle::NuMuBar)
// (206, I3Particle::NuTauBar)
// (1001, I3Particle::Brems)
// (1002, I3Particle::DeltaE)
// (1003, I3Particle::PairProd)
// (1004, I3Particle::NuclInt)
// (1005, I3Particle::MuPair)
// (1006, I3Particle::Hadrons);

inline I3Particle
to_I3Particle(const PROPOSAL::DynamicData& pp)
{
	I3Particle p;

    double x = pp.GetPosition().GetX() * I3Units::cm;
    double y = pp.GetPosition().GetY() * I3Units::cm;
    double z = pp.GetPosition().GetZ() * I3Units::cm;

    double theta = pp.GetDirection().GetTheta() * I3Units::degree;
    double phi = pp.GetDirection().GetPhi() * I3Units::degree;

    p.SetPos(x, y, z);
    p.SetThetaPhi(theta, phi);
    p.SetLength(pp.GetPropagatedDistance() * I3Units::cm);
    p.SetTime(pp.GetTime() * I3Units::second);

	p.SetType(MuonPropagator::GenerateI3Type(pp));
    p.SetLocationType(I3Particle::InIce);
    p.SetEnergy(pp.GetEnergy()*I3Units::MeV);

	return p;

	// I3Particle p;
	// particle_type_conversion_t::const_iterator it =
	//     fromRDMCTable.find(abs(pp->type));
	// if (it == fromRDMCTable.end())
	// 	log_fatal("unknown RDMC code \"%i\" cannot be converted to a I3Particle::ParticleType.", pp->type);
	// else
	// 	p.SetType(it->second);
	// p.SetLocationType(I3Particle::InIce);
	// p.SetPos(pp->x*I3Units::cm, pp->y*I3Units::cm, pp->z*I3Units::cm);
	// p.SetTime(pp->t*I3Units::s);
	// p.SetThetaPhi(pp->theta*I3Units::deg, pp->phi*I3Units::deg);
	// p.SetLength(pp->l*I3Units::cm);
	// p.SetEnergy(pp->e*I3Units::MeV);
    //
	// return p;
}

I3Particle::ParticleType MuonPropagator::GenerateI3Type(const PROPOSAL::DynamicData& secondary)
{
    PROPOSAL::DynamicData::Type type = secondary.GetTypeId();

    if (type == PROPOSAL::DynamicData::Particle)
    {
        const PROPOSAL::Particle& particle = static_cast<const PROPOSAL::Particle&>(secondary);
        PROPOSAL::ParticleDef particle_def = particle.GetParticleDef();

        I3Particle::ParticleType ptype_I3;

        bimap_ParticleType::right_const_iterator proposal_iterator = I3_PROPOSAL_ParticleType_bimap.right.find(particle_def.name);
        if (proposal_iterator == I3_PROPOSAL_ParticleType_bimap.right.end())
        {
            log_fatal("The PROPOSAL Particle '%s' can not be converted to a I3Particle", particle_def.name.c_str());
        }
        else
        {
            return proposal_iterator->second;
        }
    }
    else if(type == PROPOSAL::DynamicData::Brems) return I3Particle::Brems;
    else if(type == PROPOSAL::DynamicData::Epair) return I3Particle::PairProd;
    else if(type == PROPOSAL::DynamicData::DeltaE) return I3Particle::DeltaE;
    else if(type == PROPOSAL::DynamicData::NuclInt) return I3Particle::NuclInt;
    else
    {
        log_fatal("PROPOSAL Particle can not be converted to a I3Particle");
    }
}

/** Differential stochastic rate: d^2N/dv/dx [1/m] */
double
MuonPropagator::GetStochasticRate(double energy, double fraction, I3Particle::ParticleType type) const
{
	// propagator_->get_output()->initDefault(0, 0, GetMMCName(type), 0, 0, 0, 0, 0, 0);
    //
	// // Check kinematics
	// if (fraction <= 0 || energy*(1-fraction) <= propagator_->get_particle()->m*I3Units::MeV)
	// 	return 0.;
	// propagator_->get_particle()->setEnergy(energy/I3Units::MeV);
	// propagator_->get_cros()->get_ionization()->setEnergy();
    //
	// double contrib;
	// double rate = 0.;
	// // Separate contributions from each element for brems/epair/photonuclear interactions
	// for (int i=0; i < propagator_->get_cros()->get_medium()->get_numCompontents(); i++) {
	// 	propagator_->get_cros()->set_component(i);
	// 	if (std::isfinite(contrib = propagator_->get_cros()->get_bremsstrahlung()->get_Stochastic()->function(fraction)) && contrib > 0)
	// 		rate += contrib;
	// 	if (std::isfinite(contrib = propagator_->get_cros()->get_epairproduction()->get_Stochastic()->function(fraction)) && contrib > 0)
	// 		rate += contrib;
	// 	if (std::isfinite(contrib = propagator_->get_cros()->get_photonuclear()->get_Stochastic()->function(fraction)) && contrib > 0)
	// 		rate += contrib;
	// }
	// // Only one bulk ionization contribution
	// if (std::isfinite(contrib = propagator_->get_cros()->get_ionization()->get_Stochastic()->function(fraction)) && contrib > 0)
	// 	rate += contrib;
	// // printf("brems dN/dx: %e\n", propagator_->get_cros()->get_bremsstrahlung()->get_Stochastic()->dNdx());
	// // printf("epair dN/dx: %e\n", propagator_->get_cros()->get_epairproduction()->get_Stochastic()->dNdx());
	// // printf("photo dN/dx: %e\n", propagator_->get_cros()->get_photonuclear()->get_Stochastic()->dNdx());
	// // printf("ioniz dN/dx: %e\n", propagator_->get_cros()->get_ionization()->get_Stochastic()->dNdx());
    //
	// return rate*(I3Units::m/I3Units::cm);
}

/** total stochastic rate: dN/dx [1/m] */
double
MuonPropagator::GetTotalStochasticRate(double energy, I3Particle::ParticleType type) const
{
	// propagator_->get_output()->initDefault(0, 0, GetMMCName(type), 0, 0, 0, 0, 0, 0);
	// propagator_->get_particle()->setEnergy(energy/I3Units::MeV);
	// propagator_->get_cros()->get_ionization()->setEnergy();
    //
	// double rate = 0;
	// rate += propagator_->get_cros()->get_bremsstrahlung()->get_Stochastic()->dNdx();
	// rate += propagator_->get_cros()->get_epairproduction()->get_Stochastic()->dNdx();
	// rate += propagator_->get_cros()->get_photonuclear()->get_Stochastic()->dNdx();
	// rate += propagator_->get_cros()->get_ionization()->get_Stochastic()->dNdx();
    //
	// return rate*(I3Units::m/I3Units::cm);
}

I3Particle
MuonPropagator::propagate(const I3Particle &p, double distance, boost::shared_ptr<std::vector<I3Particle> > losses)
{
	I3Particle endpoint(p);

    double x, y, z, theta, phi = 0.0;

    PROPOSAL::Particle pp = propagator_->GetParticle();
    pp.SetParentParticleId(0);
    pp.SetParticleId(0);
    pp.SetTime(p.GetTime()/I3Units::second);

    x = p.GetPos().GetX() / I3Units::cm;
    y = p.GetPos().GetY() / I3Units::cm;
    z = p.GetPos().GetZ() / I3Units::cm;

    pp.SetPosition(PROPOSAL::Vector3D(x, y, z));
    pp.SetEnergy(p.GetEnergy()/I3Units::MeV);

    propagator_->Propagate(distance/I3Units::cm);

	endpoint.SetEnergy(pp.GetEnergy()*I3Units::MeV);

    x = pp.GetPosition().GetX() * I3Units::cm;
    y = pp.GetPosition().GetY() * I3Units::cm;
    z = pp.GetPosition().GetZ() * I3Units::cm;

    theta = pp.GetDirection().GetTheta() * I3Units::degree;
    phi = pp.GetDirection().GetPhi() * I3Units::degree;

    endpoint.SetPos(x, y, z);
    endpoint.SetThetaPhi(theta, phi);
    endpoint.SetLength(pp.GetPropagatedDistance() * I3Units::cm);
    endpoint.SetTime(pp.GetTime() * I3Units::second);

    if (losses) {
      std::vector<PROPOSAL::DynamicData*> history =
          PROPOSAL::Output::getInstance().GetSecondarys();
      for (unsigned int i = 0; i < history.size(); i++) {
        losses->push_back(to_I3Particle(*history[i]));
      }
	}

    PROPOSAL::Output::getInstance().ClearSecondaryVector();

	return endpoint;
}

void Crust::AddLayer(I3Surfaces::SurfacePtr s, boost::shared_ptr<MuonPropagator> p)
{
	boundaries_.push_back(s);
	propagators_.push_back(p);
}

I3Particle
Crust::Ingest(const I3Particle &p)
{
	I3Particle propped(p);
	double l = 0;
	for (unsigned i = 0; (propped.GetEnergy() > 0) && (i < boundaries_.size()); i++) {
		double dx = boundaries_[i]->GetIntersection(propped.GetPos(), propped.GetDir()).first;
		if (dx > 0)
			propped = (i > 0 ? propagators_[i-1] : defaultPropagator_)->propagate(propped, dx);
		// Force lengths to measure the distance back to the outermost surface
		if (i > 0)
			l += std::min(dx, propped.GetLength());
	}
	propped.SetLength(l);

	return propped;
}


}
