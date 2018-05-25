/**
 * class: I3PropagatorServicePROPOSAL
 *
 * (c) 2013 IceCube Collaboration
 */

#pragma once

// #include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>

#include "icetray/I3Logging.h"
#include "icetray/I3PointerTypedefs.h"
#include "sim-services/I3PropagatorService.h"
#include "simclasses/I3MMCTrack.h"

#include "PROPOSAL/PROPOSAL.h"

class I3Particle;

/**
 * @version $Id: I3PropagatorService.h 68823 2010-11-22 15:14:19Z dima $
 *
 * @brief this module propagates leptons through the Earth.
 *
 * @author olivas
 */

namespace PROPOSAL {

class I3PropagatorServicePROPOSAL : public I3PropagatorService
{
public:
    SET_LOGGER("I3PropagatorService");

public:

    I3PropagatorServicePROPOSAL(
        std::string configfile = "",
        I3Particle::ParticleType final_loss = I3Particle::unknown,
        double distance = 1e20
    );

    virtual ~I3PropagatorServicePROPOSAL();

    virtual std::vector<I3Particle> Propagate(I3Particle& p, DiagnosticMapPtr frame, I3FramePtr);
    virtual void SetRandomNumberGenerator(I3RandomServicePtr random);
    virtual void RegisterParticleType(I3Particle::ParticleType);

    static std::string GetDefaultConfigFile();

private:
    /** @brief Tear down and re-initialize the propagator on every call to Propagate().
     * This is used for tests that ensure that propagation results do not depend
     * on event-to-event state.
     */

    I3RandomServicePtr rng_;
    std::string config_file_;
    PropagatorService proposal_service_;

    I3Particle::ParticleType final_stochastic_loss_;
    double distance_to_propagate_;

    // default, assignment, and copy constructor declared private
    // I3PropagatorServicePROPOSAL();
    I3PropagatorServicePROPOSAL(const I3PropagatorServicePROPOSAL&);
    I3PropagatorServicePROPOSAL& operator=(const I3PropagatorServicePROPOSAL&);

    /**
     *Fills an I3MMCTrackList with position, time, and energy
     *information if the particle is a muon or a tau.
     */
    boost::shared_ptr<I3MMCTrack> GenerateMMCTrack(PROPOSAL::Particle* particle);

    boost::shared_ptr<I3MMCTrack> propagate(I3Particle& p, std::vector<I3Particle>& daughters);
};

I3_POINTER_TYPEDEFS(I3PropagatorServicePROPOSAL);
} // namespace PROPOSAL
