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
	// Kelner, Kokoulin, and Petrukhin parametrization
    int bsform = 1;
	// Abramowicz Levin Levy Maor parametrization
    int ph_fam = 3;
	// ALLM 97 (rather than 91)
    int ph_param = 2;
	// Butkevich- Mikhailov nuclear structure function
    int ph_shad = 2;

	std::ostringstream prefix;
	prefix << getenv("I3_BUILD") << "/MuonGun/resources/tables/icecube";
	//propagator_->interpolate("all", prefix.str());

    //Implement a function old param -> new param
    int new_ph_param=12;
    propagator_ = new Propagator(med,cutset,"mu",prefix.str(),molieScat,contiCorr,exactTime,lpm,1,new_ph_param,1.,1.,1.,1.,false,0);
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

typedef std::map<int, I3Particle::ParticleType> particle_type_conversion_t;

static const particle_type_conversion_t fromRDMCTable =
boost::assign::list_of<std::pair<int, I3Particle::ParticleType> >
(-100, I3Particle::unknown)
(1, I3Particle::Gamma)
(2, I3Particle::EPlus)
(3, I3Particle::EMinus)
(4, I3Particle::Nu)
(5, I3Particle::MuPlus)
(6, I3Particle::MuMinus)
(7, I3Particle::Pi0)
(8, I3Particle::PiPlus)
(9, I3Particle::PiMinus)
(11, I3Particle::KPlus)
(12, I3Particle::KMinus)
(14, I3Particle::PPlus)
(15, I3Particle::PMinus)
(33, I3Particle::TauPlus)
(34, I3Particle::TauMinus)
(41, I3Particle::Monopole)
(201, I3Particle::NuE)
(202, I3Particle::NuMu)
(203, I3Particle::NuTau)
(204, I3Particle::NuEBar)
(205, I3Particle::NuMuBar)
(206, I3Particle::NuTauBar)
(1001, I3Particle::Brems)
(1002, I3Particle::DeltaE)
(1003, I3Particle::PairProd)
(1004, I3Particle::NuclInt)
(1005, I3Particle::MuPair)
(1006, I3Particle::Hadrons);

inline I3Particle
to_I3Particle(const PROPOSALParticle *pp)
{
	I3Particle p;
	particle_type_conversion_t::const_iterator it =
	    fromRDMCTable.find(abs(pp->GetType()));
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
	propagator_->get_output()->initDefault(0, 0, GetMMCName(type), 0, 0, 0, 0, 0, 0);
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
	propagator_->get_output()->initDefault(0, 0, GetMMCName(type), 0, 0, 0, 0, 0, 0);
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
