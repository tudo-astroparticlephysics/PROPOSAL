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

#include <memory>
#include <string>
#include <vector>

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/propagation_utility/PropagationUtility.h"

namespace PROPOSAL {

class Density_distr;
class Geometry;
class Vector3D;

using Sector = std::tuple<std::shared_ptr<const Geometry>, PropagationUtility,
            std::shared_ptr<const Density_distr>>;

class Secondaries {

public:
    Secondaries(std::shared_ptr<ParticleDef> p_def, std::vector<Sector> sectors);

    // Operational functions to fill and access track
    void reserve(size_t number_secondaries);
    void clear() { track_.clear(); types_.clear(); target_hashes_.clear(); };
    void push_back(const ParticleState& point, const InteractionType& type,
                   const size_t& target_hash = 0);
    void emplace_back(const ParticleType& particle_type, const Vector3D& position,
        const Vector3D& direction, const double& energy, const double& time,
        const double& distance, const InteractionType& interaction_type,
        const size_t& target_hash = 0);
    const ParticleState& back() const {return track_.back(); }
    const ParticleState& operator[](std::size_t idx) { return track_[idx]; };

    // Track functions
    double GetELost(const Geometry& geometry) const;
    //TODO: These methods should return unique_ptr instead of shared_ptr, but pybind11 seems to have problems with them
    std::shared_ptr<ParticleState> GetEntryPoint(const Geometry& geometry) const;
    std::shared_ptr<ParticleState> GetExitPoint(const Geometry& geometry) const;
    std::shared_ptr<ParticleState> GetClosestApproachPoint(const Geometry& geometry) const;

    /*!
    * Check if particle has hit a geometry. Particle tracks ending at the border
    * of a geometry count as a hit.
    */
    bool HitGeometry(const Geometry& geometry) const;

    std::vector<ParticleState> GetTrack() const { return track_; };
    std::vector<ParticleState> GetTrack(const Geometry& geometry) const;

    ParticleState GetStateForEnergy(double energy) const;
    ParticleState GetStateForDistance(double propagated_distance) const;

    std::vector<Cartesian3D> GetTrackPositions() const;
    std::vector<Cartesian3D> GetTrackDirections() const;
    std::vector<double> GetTrackEnergies() const;
    std::vector<double> GetTrackTimes() const;
    std::vector<double> GetTrackPropagatedDistances() const;
    std::vector<InteractionType> GetTrackTypes() const { return types_; };
    std::vector<size_t> GetTargetHashes() const { return target_hashes_; };
    unsigned int GetTrackLength() const { return track_.size(); };
    std::vector<ParticleState> GetDecayProducts() const;

    // Loss functions

    std::vector<StochasticLoss> GetStochasticLosses() const;
    std::vector<StochasticLoss> GetStochasticLosses(const Geometry& geometry) const;
    std::vector<StochasticLoss> GetStochasticLosses(const InteractionType& interaction_type) const;
    std::vector<StochasticLoss> GetStochasticLosses(const std::string& interaction_type) const;

    std::vector<ContinuousLoss> GetContinuousLosses() const;
    std::vector<ContinuousLoss> GetContinuousLosses(const Geometry& geometry) const;

private:
    ParticleState RePropagateDistance(const ParticleState& init_state,
                                      const Cartesian3D& direction,
                                      double displacement) const;
    ParticleState RePropagateEnergy(const ParticleState& init_state,
                                    const Cartesian3D& direction,
                                    double energy_lost,
                                    double max_distance) const;
    Sector GetCurrentSector(const Vector3D& position,
                            const Vector3D& direction) const;

    std::vector<ParticleState> track_;
    std::vector<InteractionType> types_;
    std::vector<size_t> target_hashes_;
    std::shared_ptr<ParticleDef> primary_def_;
    std::vector<Sector> sectors_;
};

} // namespace PROPOSAL
