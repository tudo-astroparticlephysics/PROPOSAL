/** $Id: MuonPropagator.cxx 137064 2015-08-31 18:24:47Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 137064 $
 * $Date: 2015-08-31 13:24:47 -0500 (Mo, 31 Aug 2015) $
 */

#include <boost/assign.hpp>
#include <boost/foreach.hpp>

#include "MuonGun/MuonPropagator.h"

#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/PROPOSALParticle.h"
#include "PROPOSAL/Medium.h"

namespace I3MuonGun {

MuonPropagator::MuonPropagator(const std::string &medium, double ecut, double vcut, double rho)
{
	EnergyCutSettings* cutset = new EnergyCutSettings(ecut,vcut);
    Medium* med = new Medium(medium,rho);


    bool sdec      = true; // stopped muon decay
    bool exactTime = true; // exact local time
    bool molieScat = true; // Moliere scattering

	// Turn on continuous randomization if no absolute
	// energy cutoff set
    bool contiCorr = ecut < 0;

	// LPM suppression
    bool lpm = true;
    // Bremsstrahlung: Kelner, Kokoulin, and Petrukhin parametrization
    ParametrizationType::Enum brems_param = ParametrizationType::BremsKelnerKokoulinPetrukhin;
    // Photonuclear: Abramowicz Levin Levy Maor 97 with Butkevich shadowing
    ParametrizationType::Enum photo_param = ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich;

    double brems_multiplier = 1.;
    double photo_multiplier = 1.;
    double ioniz_multiplier = 1.;
    double epair_multiplier = 1.;
    bool integrate = false;
    int scattering_model = 0;

	std::ostringstream prefix;
	prefix << getenv("I3_BUILD") << "/MuonGun/resources/tables/icecube";
	//propagator_->interpolate("all", prefix.str());

    propagator_ = new Propagator(med
        , cutset
        , ParticleType::MuMinus
        , prefix.str()
        , molieScat
        , contiCorr
        , exactTime
        , lpm
        , brems_param
        , photo_param
        , brems_multiplier
        , photo_multiplier
        , ioniz_multiplier
        , epair_multiplier
        , integrate
        , scattering_model);
}

MuonPropagator::~MuonPropagator()
{
	delete propagator_;
}

void
MuonPropagator::SetSeed(int seed)
{
	MathModel::set_seed(seed);
}

ParticleType::Enum
GetPROPOSALType(I3Particle::ParticleType pt)
{
    ParticleType::Enum code;
    switch (pt) {
        case I3Particle::MuMinus:
            code = ParticleType::MuMinus;
            break;
        case I3Particle::MuPlus:
            code = ParticleType::MuPlus;
            break;
        default:
            log_fatal_stream("Unsupported particle type: " << pt);
    }
    return code;
}

std::string
MuonPropagator::GetName(const I3Particle &p)
{
    return PROPOSALParticle::GetName(GetPROPOSALType(p.GetType()));
}

typedef std::map<ParticleType::Enum, I3Particle::ParticleType> particle_type_conversion_t;

static const particle_type_conversion_t fromRDMCTable =
boost::assign::list_of<std::pair<ParticleType::Enum, I3Particle::ParticleType> >
(ParticleType::EMinus, I3Particle::EPlus)
(ParticleType::EMinus, I3Particle::EMinus)
(ParticleType::MuPlus, I3Particle::MuPlus)
(ParticleType::MuMinus, I3Particle::MuMinus)
(ParticleType::TauPlus, I3Particle::TauPlus)
(ParticleType::TauMinus, I3Particle::TauMinus)

// (-100, I3Particle::unknown)
(ParticleType::Gamma, I3Particle::Gamma)
(ParticleType::Pi0, I3Particle::Pi0)
(ParticleType::PiPlus, I3Particle::PiPlus)
(ParticleType::PiMinus, I3Particle::PiMinus)
(ParticleType::KPlus, I3Particle::KPlus)
(ParticleType::KMinus, I3Particle::KMinus)
(ParticleType::PPlus, I3Particle::PPlus)
(ParticleType::PMinus, I3Particle::PMinus)

(ParticleType::Monopole, I3Particle::Monopole)
(ParticleType::NuE, I3Particle::NuE)
(ParticleType::NuMu, I3Particle::NuMu)
(ParticleType::NuTau, I3Particle::NuTau)
(ParticleType::NuEBar, I3Particle::NuEBar)
(ParticleType::NuMuBar, I3Particle::NuMuBar)
(ParticleType::NuTauBar, I3Particle::NuTauBar)
(ParticleType::Brems, I3Particle::Brems)
(ParticleType::DeltaE, I3Particle::DeltaE)
(ParticleType::EPair, I3Particle::PairProd)
(ParticleType::NuclInt, I3Particle::NuclInt)
(ParticleType::MuPair, I3Particle::MuPair)
(ParticleType::Hadrons, I3Particle::Hadrons);

inline I3Particle
to_I3Particle(const PROPOSALParticle *pp)
{
	I3Particle p;
	particle_type_conversion_t::const_iterator it =
	    fromRDMCTable.find(pp->GetType());
	if (it == fromRDMCTable.end())
		log_fatal("unknown RDMC code \"%i\" cannot be converted to a I3Particle::ParticleType.", pp->GetType());
	else
		p.SetType(it->second);
	p.SetLocationType(I3Particle::InIce);
    p.SetPos(pp->GetX()*I3Units::cm, pp->GetY()*I3Units::cm, pp->GetZ()*I3Units::cm);
    p.SetTime(pp->GetT()*I3Units::s);
    p.SetThetaPhi(pp->GetTheta()*I3Units::deg, pp->GetPhi()*I3Units::deg);
    p.SetLength(pp->GetPropagatedDistance()*I3Units::cm);
    p.SetEnergy(pp->GetEnergy()*I3Units::MeV);

	return p;
}

/** Differential stochastic rate: d^2N/dv/dx [1/m] */
double
MuonPropagator::GetStochasticRate(double energy, double fraction, I3Particle::ParticleType type) const
{
	/*
	propagator_->get_output()->initDefault(0, 0, PROPOSALParticle::GetName(GetPROPOSALType(type)), 0, 0, 0, 0, 0, 0);
	// Check kinematics
	if (fraction <= 0 || energy*(1-fraction) <= propagator_->get_particle()->m*I3Units::MeV)
		return 0.;
	propagator_->get_particle()->setEnergy(energy/I3Units::MeV);
	propagator_->get_cros()->get_ionization()->setEnergy();
	*/

	double contrib;
	double rate = 0.;
	/*
	// Separate contributions from each element for brems/epair/photonuclear interactions
	for (int i=0; i < propagator_->get_cros()->get_medium()->get_numCompontents(); i++) {
		propagator_->get_cros()->set_component(i);
		if (std::isfinite(contrib = propagator_->get_cros()->get_bremsstrahlung()->get_Stochastic()->function(fraction)) && contrib > 0)
			rate += contrib;
		if (std::isfinite(contrib = propagator_->get_cros()->get_epairproduction()->get_Stochastic()->function(fraction)) && contrib > 0)
			rate += contrib;
		if (std::isfinite(contrib = propagator_->get_cros()->get_photonuclear()->get_Stochastic()->function(fraction)) && contrib > 0)
			rate += contrib;
	}
	// Only one bulk ionization contribution
	if (std::isfinite(contrib = propagator_->get_cros()->get_ionization()->get_Stochastic()->function(fraction)) && contrib > 0)
		rate += contrib;
	// printf("brems dN/dx: %e\n", propagator_->get_cros()->get_bremsstrahlung()->get_Stochastic()->dNdx());
	// printf("epair dN/dx: %e\n", propagator_->get_cros()->get_epairproduction()->get_Stochastic()->dNdx());
	// printf("photo dN/dx: %e\n", propagator_->get_cros()->get_photonuclear()->get_Stochastic()->dNdx());
	// printf("ioniz dN/dx: %e\n", propagator_->get_cros()->get_ionization()->get_Stochastic()->dNdx());
	*/
	return rate*(I3Units::m/I3Units::cm);
}

/** total stochastic rate: dN/dx [1/m] */
double
MuonPropagator::GetTotalStochasticRate(double energy, I3Particle::ParticleType type) const
{
	/*
	propagator_->get_output()->initDefault(0, 0, PROPOSALParticle::GetName(GetPROPOSALType(type)), 0, 0, 0, 0, 0, 0);
	propagator_->get_particle()->setEnergy(energy/I3Units::MeV);
	propagator_->get_cros()->get_ionization()->setEnergy();
	*/
	double rate = 0;
	/*
	rate += propagator_->get_cros()->get_bremsstrahlung()->get_Stochastic()->dNdx();
	rate += propagator_->get_cros()->get_epairproduction()->get_Stochastic()->dNdx();
	rate += propagator_->get_cros()->get_photonuclear()->get_Stochastic()->dNdx();
	rate += propagator_->get_cros()->get_ionization()->get_Stochastic()->dNdx();
	*/

	return rate*(I3Units::m/I3Units::cm);
}

I3Particle
MuonPropagator::propagate(const I3Particle &p, double distance, boost::shared_ptr<std::vector<I3Particle> > losses)
{
	I3Particle endpoint(p);
	/*
	// propagator_.get_output()->DEBUG=true;
	if (losses) {
		propagator_->get_output()->I3flag = true;
		propagator_->get_output()->initF2000(0, 0, GetName(p), p.GetTime()/I3Units::second,
		    p.GetPos().GetX()/I3Units::cm, p.GetPos().GetY()/I3Units::cm, p.GetPos().GetZ()/I3Units::cm,
		    p.GetDir().CalcTheta()/I3Units::deg, p.GetDir().CalcPhi()/I3Units::deg);
	} else
		propagator_->get_output()->initDefault(0, 0, GetName(p), p.GetTime()/I3Units::second,
		    p.GetPos().GetX()/I3Units::cm, p.GetPos().GetY()/I3Units::cm, p.GetPos().GetZ()/I3Units::cm,
		    p.GetDir().CalcTheta()/I3Units::deg, p.GetDir().CalcPhi()/I3Units::deg);

	PROPOSALParticle *pp = propagator_->get_particle();
	if (propagator_->propagateTo(distance/I3Units::cm, p.GetEnergy()/I3Units::MeV) > 0)
		endpoint.SetEnergy(pp->e*I3Units::MeV);
	else
		endpoint.SetEnergy(0);

	endpoint.SetPos(pp->x*I3Units::cm, pp->y*I3Units::cm, pp->z*I3Units::cm);
	endpoint.SetThetaPhi(pp->theta*I3Units::degree, pp->phi*I3Units::degree);
	endpoint.SetLength(pp->r*I3Units::cm);
	endpoint.SetTime(pp->t*I3Units::second);

	if (losses) {
		std::vector<PROPOSALParticle*> &history = propagator_->get_output()->I3hist;
		BOOST_FOREACH(PROPOSALParticle *daughter, history) {
			losses->push_back(to_I3Particle(daughter));
			delete daughter;
		}
		history.clear();
		propagator_->get_output()->I3flag = false;
	}
	*/
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
	/*
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
	*/
	return propped;
}


}
