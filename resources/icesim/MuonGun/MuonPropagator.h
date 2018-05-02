/** $Id: MuonPropagator.h 137064 2015-08-31 18:24:47Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 137064 $
 * $Date: 2015-08-31 20:24:47 +0200 (Mo, 31. Aug 2015) $
 */

#ifndef MUONGUN_MUONPROPAGATOR_H_INCLUDED
#define MUONGUN_MUONPROPAGATOR_H_INCLUDED

#include "MuonGun/I3MuonGun.h"
#include "PROPOSAL/PROPOSAL.h"
#include "dataclasses/physics/I3Particle.h"
#include "icetray/I3Units.h"
#include "phys-services/surfaces/Surface.h"

namespace I3MuonGun {

/**
 * @brief A simple muon energy-loss calculator
 *
 * This hides the nasty details of PROPOSAL (a C++ translation of MMC)
 */
class MuonPropagator
{
public:
    /**
     * @param[in] medium The name of the medium, e.g. "ice"
     * @param[in] ecut   Absolute energy above which an energy
     *                   loss is considered stochastic @f$ [MeV] @f$
     * @param[in] vcut   Proportion of the current muon energy above
     *                   which an energy loss is considered stochastic
     * @param[in] rho    Density adjustment factor for the medium
     */
    MuonPropagator(const std::string& medium, double ecut = -1, double vcut = -1, double rho = 1.0);
    ~MuonPropagator();
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
    I3Particle propagate(
        const I3Particle& p,
        double distance,
        boost::shared_ptr<std::vector<I3Particle> > losses = boost::shared_ptr<std::vector<I3Particle> >());

    /**
     * Set the (global) state of the random number generator used
     * in the implementation.
     */
    static void SetSeed(int seed);
    /**
     * Get the internal MMC name associated with a particle type
     */
    static std::string GetName(const I3Particle& p);

    double GetStochasticRate(double energy, double fraction, I3Particle::ParticleType type = I3Particle::MuMinus) const;
    double GetTotalStochasticRate(double energy, I3Particle::ParticleType type = I3Particle::MuMinus) const;

private:
    PROPOSAL::Propagator* propagator_;
    /**
     * Get the internal MMC name associated with a particle type
     */
    I3Particle::ParticleType GenerateI3Type(const PROPOSAL::DynamicData& secondary);
    I3Particle GenerateI3Particle(const PROPOSAL::DynamicData& proposal_particle);

};

/**
 * @brief A set of nested media layers
 */
class Crust
{
public:
    Crust(boost::shared_ptr<MuonPropagator> defaultPropagator)
        : defaultPropagator_(defaultPropagator){};

    /** Add an inner layer */
    void AddLayer(I3Surfaces::SurfacePtr, boost::shared_ptr<MuonPropagator>);
    /** Propagate a muon to the outer boundary of the innermost layer */
    I3Particle Ingest(const I3Particle& p);

private:
    boost::shared_ptr<MuonPropagator> defaultPropagator_;
    std::vector<I3Surfaces::SurfacePtr> boundaries_;
    std::vector<boost::shared_ptr<MuonPropagator> > propagators_;
};

} // namespace I3MuonGun

#endif // MUONGUN_MUONPROPAGATOR_H_INCLUDED
