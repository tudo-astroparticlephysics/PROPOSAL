#ifndef I3PROPAGATORSERVICE_H
#define I3PROPAGATORSERVICE_H
/**
 * class: I3PropagatorServicePROPOSAL
 *
 * (c) 2003 IceCube Collaboration
 */

#include <vector>
#include <icetray/I3Logging.h>
#include <simclasses/I3MMCTrack.h>
#include <boost/utility.hpp>
#include <icetray/I3PointerTypedefs.h>
#include <sim-services/I3PropagatorServiceBase.h>

#include "PROPOSAL/Amanda.h"
#include "PROPOSAL/PROPOSALParticle.h"
#include "PROPOSAL/Output.h"

class I3Particle; 
class PROPOSALParticle;
class Output;
class Amanda;
/**
 * @version $Id: I3PropagatorService.h 68823 2010-11-22 15:14:19Z dima $
 *
 * @brief this module propagates leptons through the Earth.
 *
 * @author olivas
 */
class I3PropagatorServicePROPOSAL : public I3PropagatorServiceBase {
 public:
 
  shared_ptr<I3MMCTrack> Propagate(I3Particle& p, std::vector<I3Particle>& daughters);
 
  SET_LOGGER("I3PropagatorService");

  /**
   * Builds an instance of this class
   * @param ctx the context with which this module's built
   */
  I3PropagatorServicePROPOSAL( const std::string&, bool debugMMC = false);

  /**
   *Destroys an instance of this class
   */
  ~I3PropagatorServicePROPOSAL(){};

 private:

  Amanda *amanda;
  Amanda *muonPropagator;
  Amanda *tauPropagator;

  double stauMass_;

  bool debugMMC_;

  std::string opts_;

  // default, assignment, and copy constructor declared private
  I3PropagatorServicePROPOSAL();
  I3PropagatorServicePROPOSAL(const I3PropagatorServicePROPOSAL&);
  I3PropagatorServicePROPOSAL& operator=(const I3PropagatorServicePROPOSAL&);

  std::string GenerateMMCName(const I3Particle&);

  /**
   *Fills an I3MMCTrackList with position, time, and energy
   *information if the particle is a muon or a tau.
   */
  boost::shared_ptr<I3MMCTrack> GenerateMMCTrack(PROPOSALParticle* particle);

  /**
   * propagate is called when flag=[1/3] to propagate given particle
   */
  shared_ptr<I3MMCTrack> propagate(I3Particle& p, std::vector<I3Particle>& daughters);

  void Fatal(const char* msg);

};

I3_POINTER_TYPEDEFS(I3PropagatorServicePROPOSAL);

#endif //I3PROPAGATORSERVICE_H
