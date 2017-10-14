
#pragma once

#include <boost/unordered_map.hpp>

#include "PROPOSAL/particle/Particle.h"

namespace PROPOSAL
{

class Propagator;

class PropagatorService
{
    public:
        typedef boost::unordered_map<ParticleDef, Propagator*> PropagatorMap;
    public:
        PropagatorService();
        virtual ~PropagatorService();

        // ----------------------------------------------------------------------------
        /// @brief Register the propagator to use
        ///
        /// The Propagators will be stored in an look up table and later called
        /// within Propagate() with the given definition of the particle
        ///
        /// @param Propagator
        // ----------------------------------------------------------------------------
        void RegisterPropagator(Propagator&);

        // ----------------------------------------------------------------------------
        /// @brief Propagate the given particle
        ///
        /// The given particle will hold the information after propagation.
        ///
        /// @param Particle
        ///
        /// @return vector of secondary data
        // ----------------------------------------------------------------------------
        std::vector<DynamicData*> Propagate(Particle&);

    private:
        PropagatorMap propagator_map_;
};

} /* PROPOSAL */

