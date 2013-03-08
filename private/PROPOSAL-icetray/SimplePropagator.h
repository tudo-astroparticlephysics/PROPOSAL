/** $Id$
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision$
 * $Date$
 */

#ifndef PROPOSAL_SIMPLEPROPAGATOR_H_INCLUDED
#define PROPOSAL_SIMPLEPROPAGATOR_H_INCLUDED

#include <boost/shared_ptr.hpp>

#include "icetray/I3Units.h"
#include "dataclasses/physics/I3Particle.h"
#include "phys-services/I3RandomService.h"

class Propagate;

namespace PROPOSAL {

/**
 * @brief A simple muon energy-loss calculator
 *
 * This hides the nasty details of PROPOSAL (a C++ translation of MMC)
 */
class SimplePropagator {
public:
	/**
	 * @param[in] medium The name of the medium, e.g. "ice"
	 * @param[in] ecut   Absolute energy above which an energy
	 *                   loss is considered stochastic @f$ [MeV] @f$
	 * @param[in] vcut   Proportion of the current muon energy above
	 *                   which an energy loss is considered stochastic
	 * @param[in] rho    Density adjustment factor for the medium
	 */
	SimplePropagator(const std::string &medium, double ecut=-1, double vcut=-1, double rho=1.0);
	~SimplePropagator();
	/**
	 * @param[in] p        Muon to propagate
	 * @param[in] distance Maximum distance to propagate
	 * @returns an I3Particle representing the muon at the end
	 *          of propagation. If the muon stopped before the
	 *          given distance, the energy will be set to zero.
	 *          The length of the output I3Particle is the
	 *          distance traveled to reach its current position
	 *          from the position given as input.
	 */
	I3Particle propagate(const I3Particle &p, double distance,
	    boost::shared_ptr<std::vector<I3Particle> > losses=boost::shared_ptr<std::vector<I3Particle> >());
	
	/**
	 * Set the (global) state of the random number generator used
	 * in the implementation.
	 */
	static void SetSeed(int seed);
	
	/**
	 * Use a specific random number generator for this instance
	 */
	void SetRandomNumberGenerator(I3RandomServicePtr rng);
	/**
	 * Get the internal MMC name associated with a particle type
	 */
	static std::string GetName(const I3Particle &p);
	
	Propagate* GetImplementation() { return propagator_; };
private:
	Propagate *propagator_;
};

}

#endif // PROPOSAL_SIMPLEPROPAGATOR_H_INCLUDED
