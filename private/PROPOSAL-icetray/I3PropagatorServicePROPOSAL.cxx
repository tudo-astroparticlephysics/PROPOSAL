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

#include <boost/foreach.hpp>
using namespace std;

class Output;

// constructor
I3PropagatorServicePROPOSAL::I3PropagatorServicePROPOSAL(const string& mmcOpts,
                                                        bool debugMMC ):
  stauMass_(NAN),
  debugMMC_(debugMMC),
  opts_(mmcOpts)
{

  //parse the opts for the stau mass
  size_t pos = mmcOpts.find("-stau"); 
  if(pos != string::npos){
    //propagating the stau and need to set the mass
    size_t begin = mmcOpts.find("=",pos) + 1;
    size_t end = mmcOpts.find(" ", begin);
    string mass = mmcOpts.substr(begin,end);
    stauMass_ = atof(mass.c_str());
  }

  
  /**
   * Configure MMC
   */


  string stderr = "/dev/null";
  if(debugMMC_)
    stderr = "";
 // if (&stderr == NULL) Fatal("Could not create the stderr string.\n");

//  env_->CallStaticBooleanMethod(Output, mid, jstr);
//  env_->DeleteLocalRef(jstr); <-----JavaStuff
 // Output::setStderr(stderr);

  amanda = new Amanda();
  amanda->setup(mmcOpts);


  log_info("PROPOSAL initialized");

}

/**
 *
 */
I3MMCTrackPtr I3PropagatorServicePROPOSAL::Propagate(I3Particle& p, vector<I3Particle>& daughters){
  // saying where we are
  log_debug("Entering I3PropagatorServicePROPOSAL::Propagate()");

  log_trace("location type = %d",p.GetLocationType());
  if(p.GetLocationType()!=I3Particle::InIce) return I3MMCTrackPtr();

  if((p.GetType()==I3Particle::NuE)||(p.GetType()==I3Particle::NuEBar)||
     (p.GetType()==I3Particle::NuMu)||(p.GetType()==I3Particle::NuMuBar)||
     (p.GetType()==I3Particle::NuTau)||(p.GetType()==I3Particle::NuTauBar)) return I3MMCTrackPtr();

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
  
  I3MMCTrackPtr mmcTrack = propagate(p, daughters);
  return mmcTrack;
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
    s<<"stau="<<stauMass_/I3Units::GeV;
    name=s.str();
  }
  //Choose the right propagator
  //if(PROP){
  //  if(p.GetType()==I3Particle::MuPlus||p.GetType()==I3Particle::MuMinus)
  //      amanda=muonPropagator;
  //  if(p.GetType()==I3Particle::TauMinus||p.GetType()==I3Particle::TauPlus)
  //      amanda=tauPropagator;
  //}
  //cout<<"amanda \t"<<amanda<<endl;
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

  //cout<<"propagate starts"<<endl;
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
  if (particle == 0) Fatal("Error calling the Particle constructor\n");

  vector<PROPOSALParticle*> aobj_l;
  aobj_l = amanda->propagate(particle);

  if (&aobj_l == NULL) {

    /**
     * If this flag flips for whatever reason this could cause mmc to fail propagation.
     */

    log_error(" Output.I3flag = %d", Output::I3flag );

    log_error("This is the properties of the particle to propagate.");
    log_error("   On the C++ side: ");
    log_error("      name = %s", mmcName.c_str() );
    log_error("      x = %f.3 cm", x_0 );
    log_error("      y = %f.3 cm ", y_0 );
    log_error("      z = %f.3 cm ", z_0 );
    log_error("      theta =%f.3 deg", theta_0 );
    log_error("      phi = %f.3 deg", phi_0 );
    log_error("      E = %f.3 MeV", e_0 );
    log_error("      t = %f.3 s", t_0 );
    log_error("   On the PROPOSAL side: ");
    log_error("      type = %d", particle->type);
    log_error("      x = %f.3 cm", particle->x );
    log_error("      y = %f.3 cm ", particle->y );
    log_error("      z = %f.3 cm ", particle->z );
    log_error("      theta =%f.3 deg",particle->theta);
    log_error("      phi = %f.3 deg", particle->phi);
    log_error("      E = %f.3 MeV", particle->e);
    log_error("      t = %f.3 s", particle->t);

    Fatal("cannot run the propagate method\n");
  }

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
    int parentID = aobj_l.at(i)->igen;
    int particleID = aobj_l.at(i)->gens;
    double x = aobj_l.at(i)->x * I3Units::cm;
    double y = aobj_l.at(i)->y * I3Units::cm;
    double z = aobj_l.at(i)->z * I3Units::cm;
    double theta = aobj_l.at(i)->theta * I3Units::deg;
    double phi = aobj_l.at(i)->phi * I3Units::deg;
    double t = aobj_l.at(i)->t * I3Units::s;
    double e = aobj_l.at(i)->e * I3Units::MeV;
    double l = aobj_l.at(i)->l * I3Units::cm;

    log_trace("MMC DEBUG SEC  \n    parent=%d pid=%d type=%d pos=(%g,%g,%g) ang=(%g,%g)  e=%g t=%g  l=%g", 
	      parentID, particleID, type, x, y, z, theta, phi, e, t, l);

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

void I3PropagatorServicePROPOSAL::Fatal(const char* msg){
  log_error("mmcOpts = %s", opts_.c_str());
  log_fatal("%s",msg);
}
