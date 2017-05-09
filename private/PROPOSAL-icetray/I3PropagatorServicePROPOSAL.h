#ifndef I3PROPAGATORSERVICEPROPOSAL_H
#define I3PROPAGATORSERVICEPROPOSAL_H
/**
 * class: I3PropagatorServicePROPOSAL
 *
 * (c) 2013 IceCube Collaboration
 */

#include <vector>
#include <boost/shared_ptr.hpp>
#include <icetray/I3Logging.h>
#include <simclasses/I3MMCTrack.h>
#include <boost/utility.hpp>
#include <icetray/I3PointerTypedefs.h>
#include <sim-services/I3PropagatorService.h>

#include "PROPOSAL/Propagator.h" //Tomasz
#include "PROPOSAL/PROPOSALParticle.h"
#include "PROPOSAL/Output.h"

class I3Particle;
class PROPOSALParticle;
class Output;
class Propagator; //Tomasz

/**
 * @version $Id: I3PropagatorService.h 68823 2010-11-22 15:14:19Z dima $
 *
 * @brief this module propagates leptons through the Earth.
 *
 * @author olivas
 */
class I3PropagatorServicePROPOSAL : public I3PropagatorService {
public:

  virtual std::vector<I3Particle> Propagate(I3Particle& p, DiagnosticMapPtr frame);

  virtual void SetRandomNumberGenerator(I3RandomServicePtr random);

  SET_LOGGER("I3PropagatorService");

  /**
   * @param mediadef[in] Path the the media definition file. If unspecified, this will
   *                     default to $I3_BUILD/PROPOSAL/resources/mediadef
   * @param tabledir[in] Path to a directory in which to store interpolation
   *                     constants for cross-section integrals. If unspecified, this will
   *                     default to $I3_BUILD/PROPOSAL/resources/tables/
   * @param cylinderRadius[in] Radius of the target volume in meters
   * @param cylinderHeight[in] Full height of the target volume in meters
   * @param type[in] Type of particle to propagate.
   * @param particleMass[in] Mass of the propagated particle in GeV. This is
   *                         only used if type is something exotic.
   * @param bs[in] Parametrization of the bremsstrahlung cross-section to use
   * @param ph[in] Family of photonuclear cross-section parametrization to use
   * @param bb[in] Specific edition of the photonuclear cross-section parametrization to use
   * @param sh[in] Nuclear shadowing parametrization to use.
   *
   * The choice of parametrizations is discussed in <a href="http://arxiv.org/abs/hep-ph/0407075">MMC paper</a>.
   */
  I3PropagatorServicePROPOSAL(std::string mediadef=""
      , std::string tabledir=""
      , double cylinderRadius=800*I3Units::m
      , double cylinderHeight=1600*I3Units::m
      , I3Particle::ParticleType type=I3Particle::MuMinus
      , double particleMass=NAN
      , ParametrizationType::Enum brems_param_ = ParametrizationType::BremsKelnerKokoulinPetrukhin
      , ParametrizationType::Enum photo_param_ = ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich
  );


  static std::string GetDefaultMediaDef();
  static std::string GetDefaultTableDir();
  virtual ~I3PropagatorServicePROPOSAL();

  void SetTearDownPerCall(bool f) { tearDownPerCall_ = f; }

 private:

  Propagator *proposal;
  double particleMass_;
  PROPOSALParticle* particle_type;

  std::string mediadef_;
  std::string tabledir_;

  double cylinderRadius_;
  double cylinderHeight_;

  ParametrizationType::Enum brems_param_;
  ParametrizationType::Enum photo_param_;

  /** @brief Tear down and re-initialize the propagator on every call to Propagate().
   * This is used for tests that ensure that propagation results do not depend
   * on event-to-event state.
   */
  bool tearDownPerCall_;
  I3RandomServicePtr rng_;

  // default, assignment, and copy constructor declared private
  // I3PropagatorServicePROPOSAL();
  I3PropagatorServicePROPOSAL(const I3PropagatorServicePROPOSAL&);
  I3PropagatorServicePROPOSAL& operator=(const I3PropagatorServicePROPOSAL&);

  ParticleType::Enum GeneratePROPOSALType(const I3Particle& p);
  I3Particle::ParticleType GenerateI3Type(ParticleType::Enum ptype_PROPOSAL);

  /**
   *Fills an I3MMCTrackList with position, time, and energy
   *information if the particle is a muon or a tau.
   */
  boost::shared_ptr<I3MMCTrack> GenerateMMCTrack(PROPOSALParticle* particle);

  boost::shared_ptr<I3MMCTrack> propagate(I3Particle& p, std::vector<I3Particle>& daughters);

};

I3_POINTER_TYPEDEFS(I3PropagatorServicePROPOSAL);

#endif //I3PROPAGATORSERVICEPROPOSAL_H
