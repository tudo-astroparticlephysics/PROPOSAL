/**
 * class: I3PropagatorServicePROPOSAL
 *
 * Mainly a copy of:
 * Version $Id: I3PropagatorService.cxx 69179 2010-12-04 01:33:36Z dima $
 *
 *
 * Date 07 Jul 2012
 *
 * (c) 2003 IceCube Collaboration
 */

// class header
#include "PROPOSAL-icetray/I3PropagatorServicePROPOSAL.h"
#include "dataclasses/physics/I3Particle.h"
#include "simclasses/I3MMCTrack.h"
#include "PROPOSAL/Amanda.h"
#include "PROPOSAL/PROPOSALParticle.h"
#include "PROPOSAL/Output.h"
#include <sstream>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>

#include <boost/foreach.hpp>
using namespace std;

class Output;

I3PropagatorServicePROPOSAL::I3PropagatorServicePROPOSAL(std::string mediadef, std::string tabledir,
    double cylinderRadius, double cylinderHeight, I3Particle::ParticleType type, double particleMass,
    BremsstrahlungParametrization bs, PhotonuclearParametrizationFamily ph, PhotonuclearParametrization bb,
    ShadowingParametrization sh) : particleMass_(particleMass)
{
	std::ostringstream mmcOpts;
	// Some options that no one should change:
	// romb: order-5 Rombert interpolation
	// raw: save interpolation tables in binary format
	// user: write "user" fields (entry/center/exit positions, energies, times for I3MMCTrack)
	// sdec: enable stopped muon decay
	// time: enable exact time calculation
	// lpm: enable Landau-Pomeranchuk-Migdal supression of EM cross-sections
	// frho: enable "smart" density factor handling
	// cont: enable continuous randomization for realistic energy loss at high vcut
	mmcOpts << "-romb=5 -raw -user -sdec -time -lpm -frho -cont";
	mmcOpts << " -bs=" << bs << " -ph=" << ph << " -bb=" << bb << " -sh=" << sh;
	mmcOpts << " -radius=" << cylinderRadius << " -height=" << cylinderHeight;
	switch (type) {
		case I3Particle::MuMinus:
		case I3Particle::MuPlus:
			break;
		case I3Particle::TauMinus:
		case I3Particle::TauPlus:
			mmcOpts << " -tau";
			break;
		case I3Particle::STauMinus:
		case I3Particle::STauPlus:
			if (!std::isfinite(particleMass))
				log_fatal("You asked for staus of non-finite mass %f", particleMass_);
			mmcOpts << " -stau=" << particleMass_;
			break;
		case I3Particle::Monopole:
			if (!std::isfinite(particleMass))
				log_fatal("You asked for monopoles of non-finite mass %f", particleMass_);
			mmcOpts << " -monopole=" << particleMass_;
			break;
		default:
			I3Particle dummy;
			dummy.SetType(type);
			log_fatal("I don't know how to propagate %s", dummy.GetTypeString().c_str());
	}
	if (mediadef.empty())
		mediadef = GetDefaultMediaDef();
	if (tabledir.empty())
		tabledir = GetDefaultTableDir();
	
	namespace fs = boost::filesystem;
	if (!fs::exists(mediadef))
		log_fatal("The mediadef file '%s' can't be read!", mediadef.c_str());
	if (!fs::is_directory(tabledir))
		log_fatal("The table directory '%s' doesn't exist!", tabledir.c_str());
	mmcOpts << " -mediadef=" << mediadef << " -tdir=" << tabledir;
	mmcOpts << " ";
	
	log_info("Amanda option string: '%s'", mmcOpts.str().c_str());
	
	amanda = new Amanda();
	amanda->setup(mmcOpts.str());
	
	mmcOpts_ = mmcOpts.str();
	tearDownPerCall_ = false;
}

std::string I3PropagatorServicePROPOSAL::GetDefaultMediaDef()
{
	const char *I3_BUILD = getenv("I3_BUILD");
	if (!I3_BUILD)
		log_fatal("$I3_BUILD is not set!");
	std::string s(I3_BUILD);
	return s + "/PROPOSAL/resources/mediadef";
}

std::string I3PropagatorServicePROPOSAL::GetDefaultTableDir()
{
	const char *I3_BUILD = getenv("I3_BUILD");
	if (!I3_BUILD)
		log_fatal("$I3_BUILD is not set!");
	std::string s(I3_BUILD);
	return s + "/PROPOSAL/resources/tables";
}

void I3PropagatorServicePROPOSAL::SetRandomNumberGenerator(I3RandomServicePtr random)
{
	rng_ = random;
	boost::function<double ()> f = boost::bind(&I3RandomService::Uniform, random, 0, 1);
	amanda->SetRandomNumberGenerator(f);
}

I3PropagatorServicePROPOSAL::~I3PropagatorServicePROPOSAL()
{
	delete amanda;
}

/**
 *
 */
std::vector<I3Particle> I3PropagatorServicePROPOSAL::Propagate(I3Particle& p, I3FramePtr frame){
  // saying where we are
  log_debug("Entering I3PropagatorServicePROPOSAL::Propagate()");
  
  vector<I3Particle> daughters;

  log_trace("location type = %d",p.GetLocationType());
  if(p.GetLocationType()!=I3Particle::InIce) return std::vector<I3Particle>();

  if((p.GetType()==I3Particle::NuE)||(p.GetType()==I3Particle::NuEBar)||
     (p.GetType()==I3Particle::NuMu)||(p.GetType()==I3Particle::NuMuBar)||
     (p.GetType()==I3Particle::NuTau)||(p.GetType()==I3Particle::NuTauBar)) return std::vector<I3Particle>();

  log_trace("particle to propagate:\n"
	    "type/energy[GeV]/posx[m]/posy[m]/posz[m]/theta[deg]/phi[deg]/length[m]\n"
	    "%d/%.2e/%.2f/%.2f/%.2f/%.2f/%.2f/%.2f",
	    p.GetType(),
	    p.GetEnergy ()/I3Units::GeV,
	    p.GetPos().GetX()/I3Units::m,
	    p.GetPos().GetY()/I3Units::m,
	    p.GetPos().GetZ()/I3Units::m,
	    p.GetZenith ()/I3Units::deg,
	    p.GetAzimuth ()/I3Units::deg,
	    p.GetLength()/I3Units::m);
  
	if (tearDownPerCall_) {
	    delete amanda;
	    amanda = new Amanda();
	    amanda->setup(mmcOpts_);
	    boost::function<double ()> f = boost::bind(&I3RandomService::Uniform, rng_, 0, 1);
	    amanda->SetRandomNumberGenerator(f);
	}
  I3MMCTrackPtr mmcTrack = propagate(p, daughters);

  if ((frame) && (mmcTrack)) {
    I3MMCTrackListConstPtr origMMCList = frame->Get<I3MMCTrackListConstPtr>("MMCTrackList");

    I3MMCTrackListPtr newMMCList;
    if (origMMCList) {
      // copy-construct a new list and delete the original one from the frame
      newMMCList = I3MMCTrackListPtr(new I3MMCTrackList(*origMMCList));
      frame->Delete("MMCTrackList");
    } else {
      newMMCList = I3MMCTrackListPtr(new I3MMCTrackList);
    }

    newMMCList->push_back(*mmcTrack);

    frame->Put("MMCTrackList", newMMCList);
  }

  return daughters;
}

string I3PropagatorServicePROPOSAL::GenerateMMCName(const I3Particle& p){

  string name;

  if(p.GetType()==I3Particle::MuMinus) name="mu-";
  else if(p.GetType()==I3Particle::MuPlus) name="mu+";
  else if(p.GetType()==I3Particle::NuMu) name="nu_mu";
  else if(p.GetType()==I3Particle::NuMuBar) name="~nu_mu";
  else if(p.GetType()==I3Particle::NuE) name="nu_e";
  else if(p.GetType()==I3Particle::NuEBar) name="~nu_e";
  else if(p.GetType()==I3Particle::NuTau) name="nu_tau";
  else if(p.GetType()==I3Particle::NuTauBar) name="~nu_tau";
  else if((p.GetType()==I3Particle::EMinus) && 
	  (p.GetShape()==I3Particle::TopShower)) name="e-";
  else if((p.GetType()==I3Particle::EPlus) && 
	  (p.GetShape()==I3Particle::TopShower)) name="e+";
  else if(p.GetType()==I3Particle::TauMinus) name="tau-";
  else if(p.GetType()==I3Particle::TauPlus) name="tau+";
  else if((p.GetType()==I3Particle::Hadrons) &&
	  (p.GetShape()==I3Particle::TopShower)) name="hadr";
  else if(p.GetType()==I3Particle::Monopole) name="monopole";
  else if(p.GetType()==I3Particle::STauPlus ||
	  p.GetType()==I3Particle::STauMinus){
    stringstream s;
    s<<"stau="<<particleMass_/I3Units::GeV;
    name=s.str();
  }

  return name;
}

I3MMCTrackPtr I3PropagatorServicePROPOSAL::GenerateMMCTrack(PROPOSALParticle* particle){
   
  //explicitly specifying the units from MMC
  double xi = particle->xi * I3Units::m;
  double yi = particle->yi * I3Units::m;
  double zi = particle->zi * I3Units::m;
  double ti = particle->ti * I3Units::second;
  double Ei = particle->Ei * I3Units::GeV;
  
  double xf = particle->xf * I3Units::m;
  double yf = particle->yf * I3Units::m;
  double zf = particle->zf * I3Units::m;
  double tf = particle->tf * I3Units::second;
  double Ef = particle->Ef * I3Units::GeV;
    
  double xc = particle->xc * I3Units::m;
  double yc = particle->yc * I3Units::m;
  double zc = particle->zc * I3Units::m;
  double tc = particle->tc * I3Units::second;
  double Ec = particle->Ec * I3Units::GeV;

  I3MMCTrackPtr mmcTrack( new I3MMCTrack);
  mmcTrack->SetEnter(xi,yi,zi,ti,Ei);
  mmcTrack->SetCenter(xc,yc,zc,tc,Ec);
  mmcTrack->SetExit(xf,yf,zf,tf,Ef);

  double Elost = particle->Elost * I3Units::GeV;
  mmcTrack->SetDepositedEnergy(Elost);
  log_debug("Elost = %f", Elost);

  if(Elost>0) return mmcTrack;
  //return null pointer if Elost <= 0
  return I3MMCTrackPtr();
}

I3MMCTrackPtr
I3PropagatorServicePROPOSAL::propagate( I3Particle& p, vector<I3Particle>& daughters){

  /**
   * Natural units of MMC is cm, deg, MeV, and s.
   * Therefore we need to convert explicitly to 
   * MMC units before passing the propagate method
   */	    
  double x_0 = p.GetPos().GetX()/I3Units::cm;     // [cm]
  double y_0 = p.GetPos().GetY()/I3Units::cm;    // [cm]
  double z_0 = p.GetPos().GetZ()/I3Units::cm;     // [cm]
  double theta_0 = p.GetDir().CalcTheta()/I3Units::deg; // [deg]
  double phi_0 = p.GetDir().CalcPhi()/I3Units::deg;   // [deg]
  double e_0 = p.GetEnergy()/I3Units::MeV;  // [MeV]
  double t_0 = p.GetTime()/I3Units::s;     // [s]
  
  string mmcName = GenerateMMCName(p);	
  log_debug("MMC name of particle to propagate: %s",mmcName.c_str());
  
  PROPOSALParticle* particle = new PROPOSALParticle(mmcName, x_0, y_0, z_0, theta_0, phi_0, e_0, t_0);
  if (particle == 0) log_fatal("Error calling the Particle constructor");

  vector<PROPOSALParticle*> aobj_l = amanda->propagate(particle);

  //get the propagated length of the particle
  double length = particle->r;
  p.SetLength( length * I3Units::cm );
  log_trace(" length = %f cm ", length );

  int nParticles =  aobj_l.size();
  log_trace("nParticles = %d", nParticles);
  
  I3MMCTrackPtr mmcTrack;
  log_trace("javaClass_ == AMANDA");
  mmcTrack = GenerateMMCTrack(particle);
  if(mmcTrack)
      mmcTrack->SetParticle( p );

  for(int i=0; i < nParticles; i++){

    
    //in mmc the particle relationships are stored
    int type = aobj_l.at(i)->type;
    double x = aobj_l.at(i)->x * I3Units::cm;
    double y = aobj_l.at(i)->y * I3Units::cm;
    double z = aobj_l.at(i)->z * I3Units::cm;
    double theta = aobj_l.at(i)->theta * I3Units::deg;
    double phi = aobj_l.at(i)->phi * I3Units::deg;
    double t = aobj_l.at(i)->t * I3Units::s;
    double e = aobj_l.at(i)->e * I3Units::MeV;
    double l = aobj_l.at(i)->l * I3Units::cm;

    log_trace("MMC DEBUG SEC  \n    type=%d pos=(%g,%g,%g) ang=(%g,%g)  e=%g t=%g  l=%g", 
	      type, x, y, z, theta, phi, e, t, l);

    //this should be a stochastic
    I3Particle new_particle;

    //Setting MC Type     
    I3Particle::ParticleType mcType(static_cast<I3Particle::ParticleType>(abs(type)));
    new_particle.SetRDMCType(mcType);
    new_particle.SetLocationType(I3Particle::InIce);
    new_particle.SetPos(x, y, z);
    new_particle.SetTime(t);
    new_particle.SetLength(l);
    new_particle.SetThetaPhi(theta,phi);
    new_particle.SetEnergy(e);

    // this is not the particle you're looking for
    // move along...and add it to the daughter list
    daughters.push_back(new_particle);

    //we're done with eobj
    delete aobj_l.at(i);
  }  

  delete particle;
  
  return mmcTrack;
}

