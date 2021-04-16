
#include "PROPOSAL/Secondaries.h"
#include "PROPOSAL/decay/DecayChannel.h"
#include "PROPOSAL/geometry/Geometry.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/density_distr/density_distr.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/Cartesian3D.h"

#include <memory>
#include <vector>

using std::get;
using namespace PROPOSAL;

Secondaries::Secondaries(std::shared_ptr<ParticleDef> p_def,
                         std::vector<Sector> sectors)
    : primary_def_(p_def)
    , sectors_(sectors)
{
}


void Secondaries::reserve(size_t number_secondaries)
{
    track_.reserve(number_secondaries);
    types_.reserve(number_secondaries);
}

void Secondaries::push_back(const ParticleState& point, const InteractionType& type)
{
    track_.push_back(point);
    types_.push_back(type);
}

void Secondaries::emplace_back(const ParticleType& particle_type, const Vector3D& position,
    const Vector3D& direction, const double& energy, const double& time,
    const double& distance, const InteractionType& interaction_type)
{
    track_.emplace_back(particle_type, position, direction, energy, time,
                        distance);
    types_.emplace_back(interaction_type);
}

std::vector<ParticleState> Secondaries::GetDecayProducts() const
{
    assert(track_.size() == types_.size());

    //TODO: Is this necessary, or do we assume that there is only one decay at the end of the vector?
    std::vector<ParticleState> decay_products;
    for (unsigned int i=0; i<track_.size(); i++) {
        if (types_[i] == InteractionType::Decay) {
            ParticleState decaying_particle = track_[i];
            double random_ch = RandomGenerator::Get().RandomDouble();
            auto products
                = primary_def_->decay_table.SelectChannel(random_ch).Decay(
                    *primary_def_, decaying_particle);
            for (auto p : products) {
                decay_products.emplace_back(p);
            }
        }
    }
    return decay_products;
}

std::vector<ParticleState> Secondaries::GetTrack(const Geometry& geometry) const
{
    std::vector<ParticleState> vec;
    for (auto i : track_) {
        if (geometry.IsInside(i.position, i.direction))
            vec.push_back(i);
    }
    return vec;
}


ParticleState Secondaries::GetStateForEnergy(double energy) const
{
    if (energy >= track_.front().energy)
        return track_.front();

    for (unsigned int i=1; i<track_.size(); i++) {
        if (track_[i].energy < energy) {
            if (types_[i] == InteractionType::ContinuousEnergyLoss) {
                return RePropagateEnergy(
                        track_[i-1], track_[i-1].energy - energy,
                        track_[i-1].propagated_distance - track_[i].propagated_distance);
            } else {
                return track_[i-1];
            }
        }
    }

    return track_.back();
}

ParticleState Secondaries::GetStateForDistance(double propagated_distance) const
{
    if (track_.front().propagated_distance >= propagated_distance)
        return track_.front();

    for (unsigned int i=1; i<track_.size(); i++) {
        if (track_[i].propagated_distance > propagated_distance) {
            return RePropagateDistance(
                    track_[i-1],
                    propagated_distance - track_[i-1].propagated_distance);
        }
    }

    return track_.back();
}

std::vector<Cartesian3D> Secondaries::GetTrackPositions() const
{
    std::vector<Cartesian3D> vec;
    for (auto i : track_)
        vec.emplace_back(i.position);
    return vec;
}

std::vector<Cartesian3D> Secondaries::GetTrackDirections() const
{
    std::vector<Cartesian3D> vec;
    for (auto i : track_)
        vec.emplace_back(i.direction);
    return vec;
}

std::vector<double> Secondaries::GetTrackEnergies() const
{
    std::vector<double> vec;
    for (auto i : track_)
        vec.emplace_back(i.energy);
    return vec;
}

std::vector<double> Secondaries::GetTrackTimes() const
{
    std::vector<double> vec;
    for (auto i : track_)
        vec.emplace_back(i.time);
    return vec;
}

std::vector<double> Secondaries::GetTrackPropagatedDistances() const
{
    std::vector<double> vec;
    for (auto i : track_)
        vec.emplace_back(i.propagated_distance);
    return vec;
}

double Secondaries::GetELost(const Geometry& geometry) const
{
    auto entry_point = GetEntryPoint(geometry);
    auto exit_point = GetExitPoint(geometry);
    return entry_point->energy - exit_point->energy;
}

std::shared_ptr<ParticleState> Secondaries::GetEntryPoint(
        const Geometry& geometry) const
{
    auto pos_0 = track_.front().position;
    auto dir_0 = track_.front().direction;
    if (geometry.IsEntering(pos_0, dir_0))
        return std::make_unique<ParticleState>(track_.front());
    if (geometry.IsInside(pos_0, dir_0))
        return nullptr; // track starts in geometry

    for (unsigned int i = 0; i < track_.size() - 1; i++) {
        auto pos_i = track_.at(i).position;
        auto pos_f = track_.at(i+1).position;
        auto dir_i = track_.at(i).direction;
        auto dist_i_f = (pos_f - pos_i).magnitude();

        auto distance = geometry.DistanceToBorder(pos_i, dir_i).first;
        if (distance <= dist_i_f && distance >= 0) {
            if ( dist_i_f - distance < GEOMETRY_PRECISION)
                return std::make_unique<ParticleState>(track_.at(i + 1));
            auto entry_point = RePropagateDistance(track_.at(i), distance);
            return std::make_unique<ParticleState>(entry_point);
        }
    }
    return nullptr; // No entry point found
}

std::shared_ptr<ParticleState> Secondaries::GetExitPoint(
        const Geometry &geometry) const
{
    auto pos_end = track_.back().position;
    auto dir_end = track_.back().direction;
    if (geometry.IsLeaving(pos_end, dir_end))
        return std::make_unique<ParticleState>(track_.back());
    if (geometry.IsInside(pos_end, dir_end))
        return nullptr; // track ends inside geometry

    for (auto i = track_.size() - 1; i > 0; i--) {
        auto pos_i = track_.at(i-1).position;
        auto pos_f = track_.at(i).position;
        auto dir_i = track_.at(i-1).direction;
        auto dist_i_f = (pos_f - pos_i).magnitude();

        auto distance = geometry.DistanceToBorder(pos_f, -dir_i).first;
        if (distance <= dist_i_f && distance >= 0) {
            if ( dist_i_f - distance < GEOMETRY_PRECISION)
                return std::make_unique<ParticleState>(track_.at(i - 1));
            auto exit_point = RePropagateDistance(track_.at(i-1),
                                                  dist_i_f - distance);
            return std::make_unique<ParticleState>(exit_point);
        }
    }

    auto pos_0 = track_.front().position;
    auto dir_0 = track_.front().direction;
    if (geometry.IsLeaving(pos_0, dir_0))
        return std::make_unique<ParticleState>(track_.front());

    return nullptr; // No exit point found
}

std::shared_ptr<ParticleState> Secondaries::GetClosestApproachPoint(const Geometry& geometry) const
{
    for (unsigned int i = 0; i < track_.size(); i++) {
        auto sec_pos = track_.at(i).position;
        auto sec_dir = track_.at(i).direction;
        if (geometry.DistanceToClosestApproach(sec_pos, sec_dir) <= 0.) {
            if(std::abs(geometry.DistanceToClosestApproach(sec_pos, sec_dir))
                        < PARTICLE_POSITION_RESOLUTION)
                return std::make_unique<ParticleState>(track_.at(i));
            if (i == 0)
                return std::make_unique<ParticleState>(track_.front());
            auto prev_pos = track_.at(i-1).position;
            auto prev_dir = track_.at(i-1).direction;
            auto displacement = geometry.DistanceToClosestApproach(prev_pos,
                                                                   prev_dir);
            auto closest_approach = RePropagateDistance(track_.at(i-1),
                                                        displacement);
            return std::make_unique<ParticleState>(closest_approach);
        }
    }
    return std::make_unique<ParticleState>(track_.back());
}

ParticleState Secondaries::RePropagateEnergy(const ParticleState& init,
                                             double energy_lost,
                                             double max_distance) const
{
    auto current_sector = GetCurrentSector(init.position,
                                           init.direction);
    auto& utility = get<Propagator::UTILITY>(current_sector);
    auto& density = get<Propagator::DENSITY_DISTR>(current_sector);

    if (utility.collection.cont_rand != nullptr) {
        Logging::Get("proposal.secondaries")->warn("Calcuating specific points "
                                                   "for sectors where continuous randomization is activated "
                                                   "can lead to inconsistent results!");
    }
    auto advance_grammage = utility.LengthContinuous(init.energy,
                                                     init.energy - energy_lost);
    auto displacement = density->Correct(init.position, init.direction,
                                         advance_grammage, max_distance);

    auto E_f = init.energy - energy_lost;
    auto new_time = init.time + utility.TimeElapsed(
            init.energy, E_f, displacement ,density->Evaluate(init.position));
    auto new_position = init.position + init.direction * displacement;
    auto new_propagated_distance = init.propagated_distance + displacement;

    return ParticleState((ParticleType)primary_def_->particle_type, new_position,
                         init.direction, E_f, new_time, new_propagated_distance);
}

ParticleState Secondaries::RePropagateDistance(const ParticleState &init,
                                               double displacement) const
{
    auto current_sector = GetCurrentSector(init.position,
                                           init.direction);
    auto& utility = get<Propagator::UTILITY>(current_sector);
    auto& density = get<Propagator::DENSITY_DISTR>(current_sector);

    if (utility.collection.cont_rand != nullptr) {
        Logging::Get("proposal.secondaries")->warn("Calcuating specific points "
                     "for sectors where continuous randomization is activated "
                     "can lead to inconsistent results!");
    }
    auto advance_grammage = density->Calculate(init.position, init.direction,
                                               displacement);
    auto E_f = utility.EnergyDistance(init.energy, advance_grammage);
    auto new_time = init.time + utility.TimeElapsed(
            init.energy, E_f,displacement, density->Evaluate(init.position));
    auto new_position = init.position + init.direction * displacement;
    auto new_propagated_distance = init.propagated_distance + displacement;
    return ParticleState((ParticleType)primary_def_->particle_type, new_position,
                         init.direction, E_f, new_time, new_propagated_distance);
}

Sector Secondaries::GetCurrentSector(const Vector3D& position,
                                     const Vector3D& direction) const
{
    //TODO: this is essentially a duplicate of Propagator::GetCurrentSector
    auto potential_sec = std::vector<Sector const*>{};
    for (auto& sector : sectors_) {
        if (get<Propagator::GEOMETRY>(sector)->IsInside(position, direction))
            potential_sec.push_back(&sector);
    }

    if (potential_sec.empty())
        throw std::logic_error(
                "Propagator: No sector defined at current particle position.");

    auto highest_sector_iter = std::max_element(
            potential_sec.begin(), potential_sec.end(), [](Sector const* a, Sector const* b) {
                return get<Propagator::GEOMETRY>(*a)->GetHierarchy()
                       < get<Propagator::GEOMETRY>(*b)->GetHierarchy();
            });

    return **highest_sector_iter;
}

std::vector<StochasticLoss> Secondaries::GetStochasticLosses() const
{
    assert(track_.size() == types_.size());

    std::vector<StochasticLoss> losses;
    for (unsigned int i=1; i<track_.size(); i++) {
        auto interaction_type = types_[i];
        if (interaction_type != InteractionType::ContinuousEnergyLoss &&
                interaction_type != InteractionType::Decay) {
            losses.emplace_back(static_cast<int>(interaction_type),
                                track_[i-1].energy - track_[i].energy,
                                track_[i].position, track_[i].direction,
                                track_[i].time, track_[i].propagated_distance,
                                track_[i-1].energy);
        }
    }
    return losses;
}

std::vector<StochasticLoss> Secondaries::GetStochasticLosses(const Geometry& geometry) const
{
    assert(track_.size() == types_.size());

    std::vector<StochasticLoss> losses;
    for (unsigned int i=1; i<track_.size(); i++) {
        if (geometry.IsInside(track_[i].position, track_[i].direction)) {
            auto interaction_type = types_[i];
            if (interaction_type != InteractionType::ContinuousEnergyLoss &&
               interaction_type != InteractionType::Decay) {
                losses.emplace_back(static_cast<int>(interaction_type),
                                    track_[i-1].energy - track_[i].energy,
                                    track_[i].position, track_[i].direction,
                                    track_[i].time, track_[i].propagated_distance,
                                    track_[i-1].energy);
            }
        }
    }
    return losses;
}

std::vector<StochasticLoss> Secondaries::GetStochasticLosses(const InteractionType& type) const
{
    assert(track_.size() == types_.size());

    std::vector<StochasticLoss> losses;
    for (unsigned int i=1; i<track_.size(); i++) {
        auto interaction_type = types_[i];
        if (interaction_type == type) {
            losses.emplace_back(static_cast<int>(interaction_type),
                                track_[i-1].energy - track_[i].energy,
                                track_[i].position, track_[i].direction,
                                track_[i].time, track_[i].propagated_distance,
                                track_[i-1].energy);
        }
    }
    return losses;
}

std::vector<StochasticLoss> Secondaries::GetStochasticLosses(
        const std::string& interaction_type) const
{
    auto it = std::find_if(Type_Interaction_Name_Map.begin(), Type_Interaction_Name_Map.end(),
                           [&interaction_type](const std::pair<InteractionType, std::string> &p) {
        return p.second == interaction_type;
    });
    if (it == Type_Interaction_Name_Map.end())
        Logging::Get("proposal.secondaries")->critical("Could not find interaction type {}", interaction_type);

    return GetStochasticLosses(it->first);
}

std::vector<ContinuousLoss> Secondaries::GetContinuousLosses() const
{
    assert(track_.size() == types_.size());

    std::vector<ContinuousLoss> losses;
    for (unsigned int i=1; i<track_.size(); i++) {
        if (types_[i] == InteractionType::ContinuousEnergyLoss) {
            losses.emplace_back( track_[i].energy - track_[i-1].energy,
                                 track_[i-1].energy, track_[i-1].position,
                                 (track_[i].position - track_[i-1].position).magnitude(),
                                 track_[i-1].direction, track_[i].direction,
                                 track_[i-1].time, track_[i].time);
        }
    }
    return losses;
}

std::vector<ContinuousLoss> Secondaries::GetContinuousLosses(const Geometry& geometry) const
{
    assert(track_.size() == types_.size());

    //TODO: At the moment, part of the continuous losses may be missing if
    // the track points are not exactly on the geometry border
    std::vector<ContinuousLoss> losses;
    for (unsigned int i=1; i<track_.size(); i++) {
        if (geometry.IsInside(track_[i].position, track_[i].direction)) {
            if (types_[i] == InteractionType::ContinuousEnergyLoss) {
                losses.emplace_back(
                        track_[i].energy - track_[i-1].energy,
                        track_[i-1].energy, track_[i-1].position,
                        (track_[i].position - track_[i-1].position).magnitude(),
                        track_[i-1].direction, track_[i].direction,
                        track_[i-1].time, track_[i].time);
            }
        }
    }
    return losses;
}
