
/******************************************************************************
 *                                                                            *
 * This file is part of the simulation tool PROPOSAL.                         *
 *                                                                            *
 * Copyright (C) 2017 TU Dortmund University, Department of Physics,          *
 *                    Chair Experimental Physics 5b                           *
 *                                                                            *
 * This software may be modified and distributed under the terms of a         *
 * modified GNU Lesser General Public Licence version 3 (LGPL),               *
 * copied verbatim in the file "LICENSE".                                     *
 *                                                                            *
 * Modifcations to the LGPL License:                                          *
 *                                                                            *
 *      1. The user shall acknowledge the use of PROPOSAL by citing the       *
 *         following reference:                                               *
 *                                                                            *
 *         J.H. Koehne et al.  Comput.Phys.Commun. 184 (2013) 2070-2090 DOI:  *
 *         10.1016/j.cpc.2013.04.001                                          *
 *                                                                            *
 *      2. The user should report any bugs/errors or improvments to the       *
 *         current maintainer of PROPOSAL or open an issue on the             *
 *         GitHub webpage                                                     *
 *                                                                            *
 *         "https://github.com/tudo-astroparticlephysics/PROPOSAL"            *
 *                                                                            *
 ******************************************************************************/


#pragma once

#include <unordered_map>

#include "PROPOSAL/particle/Particle.h"

namespace PROPOSAL {

class Propagator;

class PropagatorService
{
public:
    typedef std::unordered_map<ParticleDef, Propagator*> PropagatorMap;

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
    void RegisterPropagator(const Propagator&);

    // ----------------------------------------------------------------------------
    /// @brief Check if a propagator is registered for the given particle definition
    ///
    /// @param ParticleDef
    ///
    /// @return bool
    // ----------------------------------------------------------------------------
    bool IsRegistered(const ParticleDef&);

    // ----------------------------------------------------------------------------
    /// @brief Propagate the given particle
    ///
    /// The given particle will hold the information after propagation.
    ///
    /// @param Particle
    ///
    /// @return vector of secondary data
    // ----------------------------------------------------------------------------
    std::vector<DynamicData*> Propagate(Particle&, double distance = 1e20);

    // ----------------------------------------------------------------------------
    /// @brief Get Propagator for a given particle
    ///
    /// The Propagator to a given particle hold the information of the propagation
    /// parameters, which is useful to get access to.
    ///
    /// @param Particle
    ///
    /// @return Propagator
    // ----------------------------------------------------------------------------
    Propagator* GetPropagatorToParticleDef(const ParticleDef&);

private:
    PropagatorMap propagator_map_;
};

} // namespace PROPOSAL
