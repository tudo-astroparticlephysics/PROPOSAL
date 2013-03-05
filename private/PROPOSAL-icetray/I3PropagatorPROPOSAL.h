#ifndef I3PROPAGATORPROPOSAL_H
#define I3PROPAGATORPROPOSAL_H
/**
 * class: I3PropagatorPROPOSAL
 *
 * Version $Id: I3PropagatorPROPOSAL.h 77049 2011-06-21 19:45:18Z olivas $
 *
 * date:  9 Nov 2003
 *
 * (c) 2003 IceCube Collaboration
 */

// I3 header files
#include "dataclasses/physics/I3Particle.h"
#include "icetray/I3Module.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3EventHeader.h"
#include "simclasses/I3MMCTrack.h"
#include "PROPOSAL/Amanda.h"
#include "PROPOSAL/PROPOSALParticle.h"
#include "PROPOSAL/Output.h"

class I3Particle;
class PROPOSALParticle;
class Output;
class Amanda;

/**
 * @version $Id: I3PropagatorPROPOSAL.h 77049 2011-06-21 19:45:18Z olivas $
 *
 * @brief this module propagates leptons through the Earth.
 *
 * uses MMC version 1.3.2 (to check version, run ammc in mmc/src)
 *
 * @author dima
 */

class I3PropagatorPROPOSAL : public I3Module{
 public:
  // constructor and destructor

  SET_LOGGER("I3PropagatorPROPOSAL");

  /**
   * Builds an instance of this class
   * @param ctx the context with which this module's built
   */
  I3PropagatorPROPOSAL(const I3Context& ctx);

  /**
   *Destroys an instance of this class
   */
  ~I3PropagatorPROPOSAL();

  // transitions
  /**
   * This module takes a configuration parameter and so it must be configured
   */
  void Configure();

  /**
   * This function just calls 'Fatal' since I haven't rigged it to be
   * reconfigurable
   */
  //void Reconfigure();

  // streams
  /**
   * Dumps the event in the frame to stdout
   */
  void DAQ(I3FramePtr frame);

  void Finish();

 private:
  // parameters
  bool fPropagatorPROPOSALdebug;
  int fPropagatorPROPOSALmode;

  std::string primaryTreeName_;
  bool propagatePrimary_;
  std::string PROPOSALInfoName_;
  std::string mediadefName_;
  std::string mediadefPath_;

  /**
   * Type of exotic particle to propagate
   */
  int exoticType_;

  /**
   * Mass of the exotic particle.  This needs to be passed to PROPOSAL.
   */
  double exoticMass_;

  // add "-tau", "-e", or "-monopole[=mass in GeV]" to propagate taus, electrons, or monopoles when mode=1
  std::string fPropagatorPROPOSALopts;
  std::string fStderr;

  // default, assignment, and copy constructor declared private
  I3PropagatorPROPOSAL();
  I3PropagatorPROPOSAL(const I3PropagatorPROPOSAL&);
  I3PropagatorPROPOSAL& operator=(const I3PropagatorPROPOSAL&);

  /**
   *Fills an I3MMCTrackList with position, time, and energy
   *information if the particle is a muon or a tau.
   */
  void FillMMCTrackList(I3MMCTrackListPtr,I3Particle&,PROPOSALParticle*);

  /**
   * adds secondaries at the end of the propagation to the main track
   */
  //void eventOut(I3Particle&);
  void eventOut(I3MCTreePtr,
        I3MCTree::iterator,
        I3MMCTrackListPtr);

  std::string GetName(const I3Particle&);

  // Wrapping calls out to Java
  /**
   * set this to redirect diagnostic output from stderr to a file
   */
  void setStderr(char*);

  /**
   * call initPROPOSAL(options, flag) to initialize the PROPOSAL module:
   *       flag=1: to run tfa/Amanda
   *       flag=2: to run gen/AtmFlux as generator
   *       flag=3: to run gen/AtmFlux as propagator
   */
  void initPROPOSAL(std::string, int);

  /**
   * delete global references to PROPOSAL classes
   */
  void deletePROPOSAL();

  /**
   * propagate is called when flag=[1/3] to propagate given particle
   */
  void propagate(std::string,
         double,
         double,
         double,
         double,
         double,
         double,
         double);

  /**
   * called when done with the processed event to release resources
   */
  void endProp();

  /**
   * gets the number of particles in the returned event
   */
  int getNum();

  /**
   * called before reading particle #i
   */
  void startRead(int);



  /**
   * start/stop event (for weight accumulation)
   */
  void setStart(int);


  Amanda* amanda;
  PROPOSALParticle* particle;
  std::vector<PROPOSALParticle*> aobj;
  int flag;
};


#endif //I3PROPAGATORPROPOSAL_H
