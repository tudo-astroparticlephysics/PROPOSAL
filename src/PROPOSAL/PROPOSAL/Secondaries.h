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
    /*!
     * Creates an instance of the Secondaries class. A secondaries class serves
     * as the output of a particle propagation, and stores all information about
     * interactions during the propagation.
     *
     * For this, it stores all intermediate particle states during propagation
     * in a track object, which is a list of ParticleState objects.
     * @param p_def ParticleDef describing the physics of the propagated
     * particle
     * @param sectors List of sectors of the original propagator. This is
     * needed if particle have to be re-propagated.
     */
    Secondaries(std::shared_ptr<ParticleDef> p_def, std::vector<Sector> sectors);

    // Particle state functions

    /*!
     * Get information on the initial state of the propagated particle
     * @return ParticleState object, describing the particle at the beginning
     * of the propagation
     */
    ParticleState GetInitialState() const { return track_.front(); }

    /*!
     * Get information on the final state of the propagated particle
     * @return ParticleState object, describing the particle after the
     * propagation
     */
    ParticleState GetFinalState() const { return track_.back(); }

    /*!
     * Get information on the state of the propagated particle when it has
     * reached a specific energy.
     *
     * If the specified energy is within a continuous loss, a continuous
     * propagation step will be simulated so that the energy of the returned
     * particle state is equal to the specified energy. Note that this only
     * works consistently if continuous randomization is disabled. Otherwise, a
     * warning will be prompted and results might be inconsistent.
     *
     * If the specified energy lies within a stochastic loss (i.e. the particle
     * energy has been bigger than the specified energy before a stochastic loss
     * but smaller than the specified energy after a stochastic loss), the
     * particle state before the stochastic loss is returned.
     *
     * If the specified energy is below the final energy of the particle,
     * the final particle state will be returned.
     *
     * If the specified energy is above the initial energy of the particle,
     * the initial particle state will be returned.
     *
     * @param energy Particle energy (in MeV and in total energies) for which
     * we want to get the particle state
     * @return ParticleState object, describing the particle state at energy
     */
    ParticleState GetStateForEnergy(double energy) const;

    /*!
     * Get information on the state of the propagated particle when it has
     * reaches a specific propagation distance.
     *
     * This will calculate a continuous propagation step based on a known
     * particle state so that the propagated distance of the returned particle
     * state is equal to the specified propagated distance.
     *
     * If the specified distance is below the propagated distance of the initial
     * particle state, the initial particle state is returned.
     *
     * If the specified distance is above the propagated distance of the final
     * particle state, the final particle state is returned.
     *
     * Note that by default, the particle starts with a propagation_distance of
     * 0 unless specified otherwise in the initial_particle state.
     *
     * @param propagated_distance Propagated distance (in cm) for which we want
     * to get the particle state
     * @return ParticleState object, describing the particle state when it has
     * propagated the propagated_distance
     */
    ParticleState GetStateForDistance(double propagated_distance) const;

    //TODO: These GetXPoint methods should return unique_ptr instead of shared_ptr, but pybind11 seems to have problems with them

    /*!
     * Get the particle state when the particle has entered a specific
     * geometry.
     *
     * If the particle already started in the specific geometry or has never
     * entered the specific geometry, a nullptr is returned.
     *
     * Otherwise, a continuous step based on a known particle state will be
     * calculated so that the returned particle state describes the state when
     * the particle entered the geometry. Note that this only works reliably
     * when continuous randomization is disabled. Otherwise, a warning will be
     * prompted and results might be inconsistent.
     *
     * @param geometry Geometry for which we want to check when the propagated
     * particle has entered it
     * @return A shared pointer to a ParticleState, describing the state of the
     * particle when it has entered the specific geometry, or a nullptr.
     */
    std::shared_ptr<ParticleState> GetEntryPoint(const Geometry& geometry) const;

    /*!
     * Get the particle state when the particle has left a specific geometry.
     *
     * If the particle has never left or never crossed the geometry, a nullptr
     * is returned.
     *
     * Otherwise, a continuous step based on a known particle state will be
     * calculated so that the returned particle state describes the state when
     * the particle exited the geometry. Note that this only works reliably
     * when continuous randomization is disabled. Otherwise, a warning will be
     * prompted and results might be inconsistent.
     *
     * @param geometry Geometry for which we want to check when the propagated
     * particle has left it
     * @return A shared pointer to a ParticleState, describing the state of the
     * particle when it has left the specific geometry, or a nullptr.
     */
    std::shared_ptr<ParticleState> GetExitPoint(const Geometry& geometry) const;

    /*!
     * Get the particle state when the particle had the smallest distance to
     * the center of a specific geometry. This will not extrapolate the particle
     * before the initial or after the final state.
     *
     * This will calculate a continuous step, based on a known particle state,
     * so that the returned particle state describes the particle at is has the
     * smallest distance to the center of the specified geometry.
     * @param geometry Geometry for which we want to calculate the closest
     * approach point
     * @return A shared pointer to a ParticleState, describing the closest
     * approach point of the particle.
     */
    std::shared_ptr<ParticleState> GetClosestApproachPoint(const Geometry& geometry) const;

    /*!
     * Check if particle has hit a specific geometry. Particle tracks ending at
     * the border of a geometry count as a hit.
     * @param geometry Geometry for which we want to check whether our particle
     * has hit it
     * @return Bool - True if geometry has been hit, false otherwise.
     */
    bool HitGeometry(const Geometry& geometry) const;

    /*!
     * Calculate energy that the particle has lost inside a specific geometry.
     * @param geometry Geometry for which we want to check how much energy that
     * particle has lost inside it.
     * @return Lost energy in MeV
     */
    double GetELost(const Geometry& geometry) const;

    /*!
     * If the particle has decayed at the end of propagation, this function
     * calculated the decay products as a list of particle states. If the
     * particle did not decay during propagation, the returned list will be
     * empty.
     * @return List of ParticleStates, describing the decay products.
     */
    std::vector<ParticleState> GetDecayProducts() const;

    // Loss functions

    /*!
     * Get all stochastic losses (i.e. losses bigger than the energy cut) of the
     * particle that occurred during propagation
     * @return List of StochasticLoss objects
     */
    std::vector<StochasticLoss> GetStochasticLosses() const;

    /*!
     * Get stochastic losses (i.e. losses bigger than the energy cut) of the
     * particle that are within a specified geometry
     * @param geometry Geometry object, find stochastic losses within this geometry
     * @return List of StochasticLoss objects
     */
    std::vector<StochasticLoss> GetStochasticLosses(const Geometry& geometry) const;

    /*!
     * Get stochastic losses (i.e. losses bigger than the energy cut) of the
     * particle that correspond to a specific interaction type
     * @param interaction_type Interaction type of stochastic loss as an
     * InteractionType object
     * @return List of StochasticLoss objects
     */
    std::vector<StochasticLoss> GetStochasticLosses(const InteractionType& interaction_type) const;

    /*!
     * Get stochastic losses (i.e. losses bigger than the energy cut) of the
     * particle that correspond to a specific interaction type
     * @param interaction_type Name of interaction type of stochastic loss
     * @return List of StochasticLoss objects
     */
    std::vector<StochasticLoss> GetStochasticLosses(const std::string& interaction_type) const;

    /*!
     * Get all continuous losses (i.e. losses smaller than the energy cut) of
     * the particle that occurred during propagation
     * @return List of ContinousLoss objects
     */
    std::vector<ContinuousLoss> GetContinuousLosses() const;

    /*!
     * Get all continuous losses (i.e. losses smaller than the energy cut) of
     * the particle that occurred during propagation
     * @param geometry Geometry object, find continuous losses within this geometry
     * @return
     */
    std::vector<ContinuousLoss> GetContinuousLosses(const Geometry& geometry) const;

    // Track functions

    /*!
     * Returns particle track. This is a list of particle states that describe
     * the particle during propagation. The first element of this list describes
     * the initial particle state, the last element of this list describes the
     * particle at the end of propagation. All intermediate elements of the list
     * describe the particle during propagation. This includes particle states
     * immediately before and after stochastic losses and particle states when
     * there has been a sector transition.
     * @return List of ParticleState objects describing the states of the
     * particle during propagation.
     */
    std::vector<ParticleState> GetTrack() const { return track_; };

    /*!
     * Returns the particle track, but only the particle states of the
     * propagated particle within a specified geometry.
     * @param geometry Geometry for which we want to find the particle state
     * objects.
     * @return List of ParticleState objects, describing the intermediate states
     * of the particle during propagation within the specified geometry.
     */
    std::vector<ParticleState> GetTrack(const Geometry& geometry) const;

    /*!
     * Returns the list of all positions of the propagated particle. The first
     * element corresponds to the position of the initial particle, the last
     * element to the final position of the particle after propagation. The
     * elements in between describe the intermediate positions of the particle
     * during propagation.
     * @return List of Cartesian3D objects, describing the particle positions
     * during propagation. Note that the units of the Cartesian3D object are in
     * cm.
     */
    std::vector<Cartesian3D> GetTrackPositions() const;

    /*!
     * Returns the list of all directions of the propagated particle. The first
     * element corresponds to the direction of the initial particle, the last
     * element to the final direction of the particle after propagation. The
     * elements in between describe the intermediate directions of the particle
     * during propagation.
     * @return List of Cartesian3D objects, describing the particle directions
     * during propagation.
     */
    std::vector<Cartesian3D> GetTrackDirections() const;

    /*!
     * Returns the list of all particle energies (total energies in MeV) of the
     * propagated particle. The first element corresponds to the energy of the
     * initial particle, the last element to the final energy of the particle
     * after propagation. The elements in between describe the energies of the
     * particle during propagation.
     * @return List of doubles, describing the particle energies during
     * propagation (in total energies, in MeV)
     */
    std::vector<double> GetTrackEnergies() const;

    /*!
     * Returns the list of all particle times (times in s) of the propagated
     * particle. The first element corresponds to the time of the initial
     * particle, the last element to the time of the particle after propagation.
     * The elements in between describe the times of the particle during
     * propagation.
     * @return List of doubles, describing the particle times during propagation
     * (in s)
     */
    std::vector<double> GetTrackTimes() const;

    /*!
     * Returns the list of the propagated distances (in cm) of the propagated
     * particle. The first element corresponds to the propagated distance of the
     * initial particle (0cm if not specified otherwise in the initial state),
     * the last element to the propagated distance of the particle after
     * propagation. The elements in between describe the propagated distances
     * during propagation
     * @return List of doubles, describing the particle propagated distances
     * during propgation (in cm)
     */
    std::vector<double> GetTrackPropagatedDistances() const;

    /*!
     * Returns a list of interaction types describing the interactions of the
     * particle during propagation.
     * @return List of InteractionType objects
     */
    std::vector<InteractionType> GetTrackTypes() const { return types_; };

    /*!
     * Returns a list of hashes, describing the media and components that the
     * particle interacted with during propagation. For continuous losses,
     * 0 is returned (since the particle did not interact with a specific
     * medium).
     *
     * For ionization losses, the hash of the medium that the particle
     * interacted with is returned. We can get the medium to the corresponding
     * hash using Medium::GetMediumForHash().
     *
     * For other stochastic losses, the hash of the component that the particle
     * interacted with is returned. We can get the component to the
     * corresponding hash using Component::GetComponentForHash().
     * @return List of size_t objects, which are hashed to Medium or Component
     * objects. These are the targets that our particle interacted with during
     * propagation.
     */
    std::vector<size_t> GetTargetHashes() const { return target_hashes_; };

    /*!
     * Length of the track objects, i.e. the number of particle states that are
     * stored in this Secondaries class.
     * @return unsigned int which describes the length of the track object.
     */
    unsigned int GetTrackLength() const { return track_.size(); };

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
