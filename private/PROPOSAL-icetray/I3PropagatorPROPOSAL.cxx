/**
 * class: I3PropagatorPROPOSAL.C
 *
 * Version $Id: I3PropagatorPROPOSAL.cxx 76466 2011-06-10 17:37:38Z olivas $
 *
 * Date 07 Nov 2003
 *
 * (c) 2003 IceCube Collaboration
 */

// class header
#include "PROPOSAL-icetray/I3PropagatorPROPOSAL.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3MCTreeUtils.h"
#include "dataclasses/geometry/I3Geometry.h"
#include "icetray/I3Tray.h"
#include "icetray/I3Frame.h"
#include "simclasses/MMCWeight.h"
#include <sstream>
using namespace std;

I3_MODULE(I3PropagatorPROPOSAL);

// constructor
I3PropagatorPROPOSAL::I3PropagatorPROPOSAL(const I3Context& ctx) :
  I3Module(ctx),
  primaryTreeName_("I3MCTree"),
  propagatePrimary_(false),
  PROPOSALInfoName_("MMCTrackList"),
  exoticType_(0),
  exoticMass_(NAN)
{
  fPropagatorPROPOSALdebug=false;
  fPropagatorPROPOSALmode=1;
  AddParameter ("mode", "PROPOSAL propagation/generation mode: 1/2/3", fPropagatorPROPOSALmode);
  fPropagatorPROPOSALopts="";
  AddParameter ("opts", "PROPOSAL configuration options", fPropagatorPROPOSALopts);
  fStderr="";
  AddParameter ("rerr", "file for stderr redirection", fStderr);
  AddParameter("PrimaryTreeName","Name of Primary Tree",primaryTreeName_);
  AddParameter("PropagatePrimary","Flag to tell PROPOSAL whether to propagate primaries or not",propagatePrimary_);
  AddParameter("MMCTrackListName","Name of MMCTrack list",PROPOSALInfoName_);
  AddParameter("ShiftParticles",
           "Shifts particles to the origin for propagation and then shifts them back to the detector frame",
           bool(false) );
  AddParameter("ExoticMass","Mass of the exotic particle to propagate",exoticMass_);
  AddParameter("ExoticParticleType","Type of exotic particle to propagate",exoticType_);
  mediadefName_ = ""; // default
  AddParameter("mediadefName","name of a mediadef file. Leave it blank for default setting",mediadefName_);
  mediadefPath_= ""; // default
  AddParameter("mediadefPath","Path to mediadef file. Set it if you define mediadefName",mediadefPath_);
  // Specifying the outbox
  AddOutBox("OutBox");
}

// destructor
I3PropagatorPROPOSAL::~I3PropagatorPROPOSAL()
{
}

void I3PropagatorPROPOSAL::Finish()
{
  deletePROPOSAL();
}

// configure
void I3PropagatorPROPOSAL::Configure(){
  // saying where we are
  log_debug("Configuring the I3PropagatorPROPOSAL Module");

  // retrieving the parameters
  GetParameter ("mode", fPropagatorPROPOSALmode);
  log_info("PROPOSAL propagation/generation mode set to %d", fPropagatorPROPOSALmode);

  GetParameter ("rerr", fStderr);
  char rerr[1000]="\0";
  strcat(rerr, fStderr.c_str());
  if(rerr[0]!='\0') setStderr(rerr);

  if(fPropagatorPROPOSALmode<0){
    fPropagatorPROPOSALmode=-fPropagatorPROPOSALmode;
    fPropagatorPROPOSALdebug=true;
  }
  else fPropagatorPROPOSALdebug=false;
  GetParameter ("opts", fPropagatorPROPOSALopts);

  GetParameter("ExoticMass",exoticMass_);
  GetParameter("ExoticParticleType",exoticType_);

  GetParameter("mediadefName",mediadefName_);
  GetParameter("mediadefPath",mediadefPath_);

  const string I3_BUILD(getenv("I3_BUILD"));
  string options;
  options = "-romb=5 -raw -user -sdec -time -lpm -bs=1 -ph=3 -bb=2 -sh=2 ";
  options += "-frho -tdir="+I3_BUILD+"/PROPOSAL/resources/tables ";

  if(fPropagatorPROPOSALmode==3) options+="-prop ";

  options += fPropagatorPROPOSALopts;

  if(exoticType_ == static_cast<int>(I3Particle::Monopole) ||
     exoticType_ == static_cast<int>(I3Particle::STauPlus) ||
     exoticType_ == static_cast<int>(I3Particle::STauMinus) ){
    stringstream s;
    if(exoticType_ == static_cast<int>(I3Particle::Monopole) ){
      s<<" -monopole="<<exoticMass_/I3Units::GeV;
    }
    if(exoticType_ == static_cast<int>(I3Particle::STauPlus) ||
       exoticType_ == static_cast<int>(I3Particle::STauMinus) ){
      s<<" -stau="<<exoticMass_/I3Units::GeV;
    }
    options += s.str();
    options += " -mediadef="+I3_BUILD+"/PROPOSAL/resources/mediadef-exotics ";
  }else{
    //don't track constant energy loss for very massive particles
    //but we do want to for everything else.
    options += " -cont ";
    if (mediadefName_ == "") {
       // default
       options += " -mediadef="+I3_BUILD+"/PROPOSAL/resources/mediadef ";
    } else {
       // user-defined media file. Put it under I3_BUILD
       if (mediadefPath_ == "") mediadefPath_ = I3_BUILD;
       options += " -mediadef="+mediadefPath_ +"/"+mediadefName_+" ";
    }
  }

  char opts[1000]="\0";
  strcat(opts, options.c_str());
  initPROPOSAL(opts, fPropagatorPROPOSALmode);

  GetParameter("PrimaryTreeName",primaryTreeName_);
  GetParameter("PropagatePrimary",propagatePrimary_);
  GetParameter("MMCTrackListName",PROPOSALInfoName_);
  bool deprecated(false);
  GetParameter("ShiftParticles",deprecated);
  if(deprecated) log_error("I3PropagatorPROPOSAL no longer shifts partciles.");

}

/**
 * gets muons/taus in mode 1 or all leptons in mode 3 and calls the propagate method,
 * filling secondaries as necessary.  Calls createNext method in mode 2 to create and
 * propagate leptons, and adds them and all secondaries to the frame.
 */
void I3PropagatorPROPOSAL::DAQ(I3FramePtr frame){
  // saying where we are
  log_debug("Entering I3PropagatorPROPOSAL::DAQ()");
  I3MMCTrackListPtr mmcTrackList(new I3MMCTrackList);
  
    //The primary tree
    I3MCTreeConstPtr primtree = I3MCTreeUtils::Get(frame,primaryTreeName_);

    if((!primtree) || (primtree->size()==0)) {
      log_debug("In propagate mode but I have nothing to propagate.");
    }else{
      log_debug("starting to process primtree '%s'",primaryTreeName_.c_str());

      I3MCTreePtr proptree(new I3MCTree(*primtree));

      if(fPropagatorPROPOSALmode==1) setStart(1);

      for(I3MCTree::post_order_iterator iter = proptree->begin_post(); iter!=proptree->end_post(); ++iter){
    //Skip primaries if that's what the user wants
    log_trace("depth = %d propagatePrimary_ = %d",proptree->depth(iter),propagatePrimary_);
    if((proptree->depth(iter) == 0) && (!propagatePrimary_)) continue;

    log_trace("location type = %d",iter->GetLocationType());
    if(iter->GetLocationType()!=I3Particle::InIce) continue;

        if((iter->GetType()==I3Particle::NuE)||(iter->GetType()==I3Particle::NuEBar)||
       (iter->GetType()==I3Particle::NuMu)||(iter->GetType()==I3Particle::NuMuBar)||
       (iter->GetType()==I3Particle::NuTau)||(iter->GetType()==I3Particle::NuTauBar)) continue;

    log_debug("particle to propagate:\n"
          "type/energy[GeV]/posx[m]/posy[m]/posz[m]/theta[deg]/phi[deg]/length[m]\n"
          "%d/%.2e/%.2f/%.2f/%.2f/%.2f/%.2f/%.2f",
          iter->GetType(),
          iter->GetEnergy ()/I3Units::GeV,
          iter->GetPos().GetX()/I3Units::m,
          iter->GetPos().GetY()/I3Units::m,
          iter->GetPos().GetZ()/I3Units::m,
          iter->GetZenith ()/I3Units::deg,
          iter->GetAzimuth ()/I3Units::deg,
          iter->GetLength()/I3Units::m);

      string name = GetName(*iter);
      log_trace("Name of particle to propagate: %s",name.c_str());
      if(name.size()){

        /**
         * Natural units of PROPOSAL is cm, deg, MeV, and s.
         * Therefore we need to convert explicitly to
         * PROPOSAL units before passing the propagate method
         */
        propagate(name.c_str(),
              iter->GetPos().GetX()/I3Units::cm,
              iter->GetPos().GetY()/I3Units::cm,
              iter->GetPos().GetZ()/I3Units::cm,
              iter->GetDir().CalcTheta()/I3Units::deg,
              iter->GetDir().CalcPhi()/I3Units::deg,
              iter->GetEnergy()/I3Units::MeV,
              iter->GetTime()/I3Units::s);

        iter->SetLength( particle->r * I3Units::cm );

        /**
         * eventOut fills the tree with the secondaries
         * We need to get the iterator in the new tree so that eventOut
         * can add children correctly.  This is faster than searching
         * for the parent by particle ID for every child.
         */
        eventOut(proptree,iter,mmcTrackList);
        if(flag==1) FillMMCTrackList(mmcTrackList,*iter,particle);
        endProp();
      }
      }// end primary iter
      if(fPropagatorPROPOSALmode==1 && amanda->rw!=0){
    setStart(0);

    log_debug("EVENT WEIGHT = %f %f", amanda->fw, amanda->hw);

    //make an MMCWeight object to hold the weight and put it in the frame...
    MMCWeightPtr weight = MMCWeightPtr(new MMCWeight());
    weight->weight = amanda->fw;
    weight->distToModIntPoint = amanda->hw;

    log_debug("Putting the PROPOSAL weight in the frame...");
    frame->Put("MMCWeight", weight);
      }
      frame->Delete(primaryTreeName_);

      frame->Put(primaryTreeName_, proptree);
    }

  frame->Put(PROPOSALInfoName_,mmcTrackList);

  // sending out the event
  log_debug("Pushing the frame");
  PushFrame(frame,"OutBox");
}

string I3PropagatorPROPOSAL::GetName(const I3Particle& p){

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
    s<<"stau="<<exoticMass_/I3Units::GeV;
    name=s.str();
  }

  return name;
}

/**
 * gets parameters of the propagated lepton, and adds them together with all secondaries to the lepton track
 */
void
I3PropagatorPROPOSAL::eventOut(I3MCTreePtr tree,
              I3MCTree::iterator iter,
              I3MMCTrackListPtr mmctracks){

  string name="mu-";
  int type;
  double x=0;     // [cm]
  double y=0;     // [cm]
  double z=0;     // [cm]
  double theta=0; // [deg]
  double phi=0;   // [deg]
  double e=1000;  // [MeV]
  double l=0;     // [cm]
  double t=0;     // [s]

  for(int i=0; i < aobj.size(); i++){

    name = aobj.at(i)->name;
    type = aobj.at(i)->type;
    x = aobj.at(i)->x * I3Units::cm;
    y = aobj.at(i)->y * I3Units::cm;
    z = aobj.at(i)->z * I3Units::cm;
    theta = aobj.at(i)->theta * I3Units::deg;
    phi = aobj.at(i)->phi * I3Units::deg;
    t = aobj.at(i)->t * I3Units::s;
    e = aobj.at(i)->e * I3Units::MeV;
    l = aobj.at(i)->r * I3Units::cm;

    // what to do for undefined l? (PROPOSAL returns -1 when particles exit the Earth)
    if(fPropagatorPROPOSALdebug) log_debug("PROPOSAL DEBUG SEC   %s  %g %g %g  %g %g  %g %g  %g",
                                           name.c_str(), x, y, z, theta, phi, e, t, l);

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
    tree->append_child(iter,new_particle);

    FillMMCTrackList(mmctracks,new_particle,aobj.at(i));
  }
}

void I3PropagatorPROPOSAL::FillMMCTrackList(I3MMCTrackListPtr t,I3Particle& p,PROPOSALParticle* par){
  if(p.GetType()==I3Particle::MuPlus ||
     p.GetType()==I3Particle::MuMinus ||
     p.GetType()==I3Particle::STauPlus ||
     p.GetType()==I3Particle::STauMinus ||
     p.GetType()==I3Particle::TauPlus ||
     p.GetType()==I3Particle::TauMinus ){

    log_debug("Creating an I3MMCInfo object for this particle.");
    I3MMCTrack mmcTrack;
    mmcTrack.SetParticle(p);

    //explicitly specifying the units from PROPOSAL
    double xi = par->xi * I3Units::m;
    double yi = par->yi * I3Units::m;
    double zi = par->zi * I3Units::m;
    double ti = par->ti * I3Units::second;
    double Ei = par->Ei * I3Units::GeV;

    double xf = par->xf * I3Units::m;
    double yf = par->yf * I3Units::m;
    double zf = par->zf * I3Units::m;
    double tf = par->tf * I3Units::second;
    double Ef = par->Ef * I3Units::GeV;

    double xc = par->xc * I3Units::m;
    double yc = par->yc * I3Units::m;
    double zc = par->zc * I3Units::m;
    double tc = par->tc * I3Units::second;
    double Ec = par->Ec * I3Units::GeV;

    mmcTrack.SetEnter(xi,yi,zi,ti,Ei);
    mmcTrack.SetCenter(xc,yc,zc,tc,Ec);
    mmcTrack.SetExit(xf,yf,zf,tf,Ef);

    double Elost = par->Elost * I3Units::GeV;
    mmcTrack.SetDepositedEnergy(Elost);

    if(Elost>0) t->push_back(mmcTrack);
  }
}

void
I3PropagatorPROPOSAL::setStderr(char *filename){}

void
I3PropagatorPROPOSAL::initPROPOSAL(string options, int iflag){

  flag=iflag;

  amanda = new Amanda();
  amanda->setup(options);

  log_info("PROPOSAL initialized");
}


void
I3PropagatorPROPOSAL::deletePROPOSAL(){
 // env->DeleteGlobalRef(amanda);
 // log_info("PROPOSAL deleted");
}


void
I3PropagatorPROPOSAL::propagate(string name,
               double x,
               double y,
               double z,
               double theta,
               double phi,
               double e,
               double t)
{


  particle = new PROPOSALParticle(name, x, y, z, theta, phi, e, t);

  // propagate the particle

  aobj = amanda->propagate(particle);

}

void
I3PropagatorPROPOSAL::endProp(){

}



void
I3PropagatorPROPOSAL::setStart(int f){
  if(amanda->rw!=0) amanda->dw = f;
}
