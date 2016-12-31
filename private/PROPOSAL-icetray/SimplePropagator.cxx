/** $Id: SimplePropagator.cxx 150208 2016-09-21 11:32:47Z tfuchs $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 150208 $
 * $Date: 2016-09-21 13:32:47 +0200 (Mi, 21. Sep 2016) $
 */

#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/bind.hpp>

#include "PROPOSAL-icetray/SimplePropagator.h"

#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Photonuclear.h"
#include "PROPOSAL/MathModel.h"

#include "PROPOSAL/Medium.h"

namespace PROPOSAL {

SimplePropagator::SimplePropagator(const std::string &medium, double ecut, double vcut, double rho)
//--- Tomasz
    //: propagator_(new Propagate(medium, ecut, vcut, "mu", rho))
{
    EnergyCutSettings* cutset = new EnergyCutSettings(ecut,vcut);
    Medium* med = new Medium(medium,rho);

    //Tomasz
//	propagator_->sdec      = true; // stopped muon decay
//	propagator_->exactTime = true; // exact local time
//	propagator_->molieScat = true; // Moliere scattering
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

SimplePropagator::~SimplePropagator()
{
	delete propagator_;
}

void
SimplePropagator::SetSeed(int seed)
{
	MathModel::set_seed(seed);
}

void
SimplePropagator::SetRandomNumberGenerator(I3RandomServicePtr rng)
{
	boost::function<double ()> f = boost::bind(&I3RandomService::Uniform, rng, 0, 1);
	propagator_->SetRandomNumberGenerator(f);
}

inline std::string
GetMMCName(I3Particle::ParticleType pt)
{
	std::string name;
	
	switch (pt) {
		case I3Particle::MuMinus:
			name="mu-";
			break;
		case I3Particle::MuPlus:
			name="mu+";
			break;
		default:
			break;
	}
	
	return name;
}

std::string
SimplePropagator::GetName(const I3Particle &p)
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
	I3Particle::ParticleType type;
	particle_type_conversion_t::const_iterator it =
	    fromRDMCTable.find(abs(pp->GetType())); //Tomasz
	if (it == fromRDMCTable.end()) {
		log_fatal("unknown RDMC code \"%i\" cannot be converted to a I3Particle::ParticleType. It will appear as \"unknown\".", pp->GetType());//Tomasz
		type = I3Particle::unknown;
	} else {
		type = it->second;
	}
	p.SetType(type);
	//Tomasz
    p.SetLocationType(I3Particle::InIce);
    p.SetPos(pp->GetX()*I3Units::cm, pp->GetY()*I3Units::cm, pp->GetZ()*I3Units::cm);
    p.SetTime(pp->GetT()*I3Units::s);
    p.SetThetaPhi(pp->GetTheta()*I3Units::deg, pp->GetPhi()*I3Units::deg);
    p.SetLength(pp->GetPropagatedDistance()*I3Units::cm);
    p.SetEnergy(pp->GetEnergy()*I3Units::MeV);
	//Tomasz End
	return p;
}

I3Particle
SimplePropagator::propagate(const I3Particle &p, double distance, boost::shared_ptr<std::vector<I3Particle> > losses)
{
	I3Particle endpoint(p);
	
	/*
	if (losses) {
		propagator_->get_output()->I3flag = true;
		propagator_->get_output()->initF2000(0, 0, GetName(p), p.GetTime()/I3Units::second,
		    p.GetPos().GetX()/I3Units::cm, p.GetPos().GetY()/I3Units::cm, p.GetPos().GetZ()/I3Units::cm,
		    p.GetDir().CalcTheta()/I3Units::deg, p.GetDir().CalcPhi()/I3Units::deg);
	} else
		propagator_->get_output()->initDefault(0, 0, GetName(p), p.GetTime()/I3Units::second,
		    p.GetPos().GetX()/I3Units::cm, p.GetPos().GetY()/I3Units::cm, p.GetPos().GetZ()/I3Units::cm,
		    p.GetDir().CalcTheta()/I3Units::deg, p.GetDir().CalcPhi()/I3Units::deg);
	
	*/
	
	    PROPOSALParticle *pp = propagator_->GetParticle();
    pp->SetParentParticleId(0);
    pp->SetParticleId(0);
    pp->SetT(p.GetTime()/I3Units::second);
    pp->SetX(p.GetPos().GetX()/I3Units::cm);
    pp->SetY(p.GetPos().GetY()/I3Units::cm);
    pp->SetZ(p.GetPos().GetZ()/I3Units::cm);
    pp->SetEnergy(p.GetEnergy()/I3Units::MeV);



    
    if (propagator_->Propagate(distance/I3Units::cm))
	    endpoint.SetEnergy(pp->GetEnergy()*I3Units::MeV);
	else
	    endpoint.SetEnergy(0);

	endpoint.SetPos(pp->GetX()*I3Units::cm, pp->GetY()*I3Units::cm, pp->GetZ()*I3Units::cm);
	endpoint.SetThetaPhi(pp->GetTheta()*I3Units::degree, pp->GetPhi()*I3Units::degree);
	endpoint.SetLength(pp->GetPropagatedDistance()*I3Units::cm);
	endpoint.SetTime(pp->GetT()*I3Units::second);
	
	//Tomasz
	if (losses) {
	    //std::vector<PROPOSALParticle*> &history = propagator_->get_output()->I3hist;
	    std::vector<PROPOSALParticle*> history = Output::getInstance().GetSecondarys();
	    for(unsigned int i = 0; i<history.size();i++)
	    {
		losses->push_back(to_I3Particle(history.at(i)));
	    }
	}
	Output::getInstance().ClearSecondaryVector();
	//Tomasz New End
	
	//Tomasz
	/*
	if (losses) {
		std::vector<PROPOSALParticle*> &history = propagator_->get_output()->I3hist;
		BOOST_FOREACH(PROPOSALParticle *pp, history) {
			losses->push_back(to_I3Particle(pp));
			delete pp;
		}
		history.clear();
		propagator_->get_output()->I3flag = false;
	}
	*/ //Tomasz Old End
	return endpoint;
}

}
