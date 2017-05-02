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
#include <boost/bimap.hpp>

#include "PROPOSAL-icetray/SimplePropagator.h"

#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/Photonuclear.h"
#include "PROPOSAL/MathModel.h"

#include "PROPOSAL/Medium.h"

namespace PROPOSAL {


ParticleType::Enum GetPROPOSALType(I3Particle::ParticleType pt)
{
    ParticleType::Enum particle_type;
    switch (pt) {
        case I3Particle::MuMinus:
            particle_type = ParticleType::MuMinus;
            break;
        case I3Particle::MuPlus:
            particle_type = ParticleType::MuPlus;
            break;
        case I3Particle::TauMinus:
            particle_type = ParticleType::TauMinus;
            break;
        case I3Particle::TauPlus:
            particle_type = ParticleType::TauPlus;
            break;
        default:
            log_fatal_stream("Unsupported particle type: " << pt);
    }
    return particle_type;
}


SimplePropagator::SimplePropagator(const std::string &medium, I3Particle::ParticleType pt, double ecut, double vcut, double rho)
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
    // TODO: Implement a function old param -> new param or just do it in enums

   double brems_multiplier = 1.;
   double photo_multiplier = 1.;
   double ioniz_multiplier = 1.;
   double epair_multiplier = 1.;
   bool integrate = false;
   int scattering_model = 0;

    std::ostringstream prefix;
    prefix << getenv("I3_BUILD") << "/MuonGun/resources/tables/icecube";
    //propagator_->interpolate("all", prefix.str());

    propagator_ = new Propagator(med,cutset
    	, GetPROPOSALType(pt)
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
    propagator_->SetStopping_decay(sdec);
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

std::string
SimplePropagator::GetName(const I3Particle &p)
{
    return PROPOSALParticle::GetName(GeneratePROPOSALType(p));
}

typedef boost::bimap<I3Particle::ParticleType, ParticleType::Enum> bimap_ParticleType;
static const bimap_ParticleType I3_PROPOSAL_ParticleType_bimap = boost::assign::list_of<bimap_ParticleType::relation>
    (I3Particle::MuMinus,   ParticleType::MuMinus)
    (I3Particle::MuPlus,    ParticleType::MuPlus)
    (I3Particle::TauMinus,  ParticleType::TauMinus)
    (I3Particle::TauPlus,   ParticleType::TauPlus)
    (I3Particle::EMinus,    ParticleType::EMinus)
    (I3Particle::EPlus,     ParticleType::EPlus)
    (I3Particle::NuMu,      ParticleType::NuMu)
    (I3Particle::NuMuBar,   ParticleType::NuMuBar)
    (I3Particle::NuE,       ParticleType::NuE)
    (I3Particle::NuEBar,    ParticleType::NuEBar)
    (I3Particle::NuTau,     ParticleType::NuTau)
    (I3Particle::NuTauBar,  ParticleType::NuTauBar)
    (I3Particle::Brems,     ParticleType::Brems)
    (I3Particle::DeltaE,    ParticleType::DeltaE)
    (I3Particle::PairProd,  ParticleType::EPair)
    (I3Particle::NuclInt,   ParticleType::NuclInt)
    (I3Particle::MuPair,    ParticleType::MuPair)
    (I3Particle::Hadrons,   ParticleType::Hadrons)
    (I3Particle::Monopole,  ParticleType::Monopole)
    (I3Particle::STauMinus, ParticleType::STauMinus)
    (I3Particle::STauPlus,  ParticleType::STauPlus);


// ------------------------------------------------------------------------- //
ParticleType::Enum SimplePropagator::GeneratePROPOSALType(const I3Particle& p)
{
    I3Particle::ParticleType ptype_I3 = p.GetType();
    ParticleType::Enum ptype_PROPOSAL;

    bimap_ParticleType::left_const_iterator i3_iterator = I3_PROPOSAL_ParticleType_bimap.left.find(ptype_I3);
    if (i3_iterator == I3_PROPOSAL_ParticleType_bimap.left.end())
    {
        log_fatal("The I3Particle '%s' with type '%i' can not be converted to a PROPOSALParticle"
            , p.GetTypeString().c_str(), ptype_I3);
    }

    ptype_PROPOSAL = I3_PROPOSAL_ParticleType_bimap.left.find(ptype_I3) -> second;

    return ptype_PROPOSAL;
}

I3Particle::ParticleType SimplePropagator::GenerateI3Type(ParticleType::Enum ptype_PROPOSAL)
{
    I3Particle::ParticleType ptype_I3;

    bimap_ParticleType::right_const_iterator proposal_iterator = I3_PROPOSAL_ParticleType_bimap.right.find(ptype_PROPOSAL);
    if (proposal_iterator == I3_PROPOSAL_ParticleType_bimap.right.end())
    {
        log_fatal("The PROPOSALParticle '%s' with type '%i' can not be converted to a I3Particle"
            , PROPOSALParticle::GetName(ptype_PROPOSAL).c_str(), ptype_PROPOSAL);
    }

    ptype_I3 = I3_PROPOSAL_ParticleType_bimap.right.find(ptype_PROPOSAL) -> second;

    return ptype_I3;
}

inline I3Particle
to_I3Particle(const PROPOSALParticle *pp)
{
	I3Particle p;

	p.SetType(GenerateI3Type(pp->GetType()));
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
