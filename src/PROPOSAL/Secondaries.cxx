
#include "PROPOSAL/Secondaries.h"
#include "PROPOSAL/decay/DecayChannel.h"
#include "PROPOSAL/geometry/Geometry.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/Logging.h"

#include <memory>
#include <vector>

using namespace PROPOSAL;

Secondaries::Secondaries(std::shared_ptr<ParticleDef> p_def,
                         std::vector<Sector> sectors)
    : primary_def_(p_def)
    , sectors_(sectors)
{
}


void Secondaries::reserve(size_t number_secondaries)
{
    secondaries_.reserve(number_secondaries);
}

void Secondaries::push_back(const DynamicData& continuous_loss)
{
    secondaries_.push_back(continuous_loss);
}

void Secondaries::emplace_back(const ParticleType& type, const Vector3D& position,
    const Vector3D& direction, const double& energy, const double& time,
    const double& distance)
{
    secondaries_.emplace_back(type, position, direction, energy, time, distance);
}

// void Secondaries::push_back(const Particle& particle, const int&
// interaction_type, const double& energy_loss)
// {
//     DynamicData data(interaction_type);

//     data.SetEnergy(energy_loss);
//     data.SetPosition(particle.GetPosition());
//     data.SetDirection(particle.GetDirection());
//     data.SetTime(particle.GetTime());
//     data.SetParentParticleEnergy(particle.GetEnergy());
//     data.SetPropagatedDistance(particle.GetPropagatedDistance());

//     secondaries_.push_back(data);
// }

void Secondaries::append(Secondaries& secondaries)
{
    secondaries_.insert(secondaries_.end(), secondaries.secondaries_.begin(),
        secondaries.secondaries_.end());
}

Secondaries Secondaries::Query(const int& interaction_type) const
{
    Secondaries sec(primary_def_, sectors_);
    for (auto i : secondaries_) {
        if (interaction_type == i.type)
            sec.push_back(i);
    }
    return sec;
}

Secondaries Secondaries::Query(const std::string& interaction_type) const
{
    Secondaries sec(primary_def_, sectors_);
    for (auto i : secondaries_) {
        if (interaction_type == i.GetName())
            sec.push_back(i);
    }
    return sec;
}

Secondaries Secondaries::Query(const Geometry& geometry) const
{
    Secondaries sec(primary_def_, sectors_);
    for (auto i : secondaries_) {
        if (geometry.IsInside(i.position, i.direction))
            sec.push_back(i);
    }
    return sec;
}

void Secondaries::DoDecay()
{
    for (auto it = secondaries_.begin(); it != secondaries_.end();) {
        if (it->type == static_cast<int>(InteractionType::Decay)) {
            DynamicData decaying_particle(static_cast<ParticleType>(primary_def_->particle_type),
                it->position, it->direction, it->energy, it->time,
                it->propagated_distance);
            double random_ch = RandomGenerator::Get().RandomDouble();
            auto products
                = primary_def_->decay_table.SelectChannel(random_ch).Decay(
                    *primary_def_, decaying_particle);
            it = secondaries_.erase(it); // delete old decay
            for (auto p : products) {
                // and insert decayparticles inplace of old decay
                it = secondaries_.insert(it, std::move(p));
                it++;
            }
        } else {
            it++;
        }
    }
}

std::vector<Vector3D> Secondaries::GetPosition() const
{
    std::vector<Vector3D> vec;
    for (auto i : secondaries_)
        vec.emplace_back(i.position);
    return vec;
}

std::vector<Vector3D> Secondaries::GetDirection() const
{
    std::vector<Vector3D> vec;
    for (auto i : secondaries_)
        vec.emplace_back(i.direction);
    return vec;
}

std::vector<double> Secondaries::GetEnergy() const
{
    std::vector<double> vec;
    for (auto i : secondaries_)
        vec.emplace_back(i.energy);
    return vec;
}



std::vector<double> Secondaries::GetTime() const
{
    std::vector<double> vec;
    for (auto i : secondaries_)
        vec.emplace_back(i.time);
    return vec;
}

std::vector<double> Secondaries::GetPropagatedDistance() const
{
    std::vector<double> vec;
    for (auto i : secondaries_)
        vec.emplace_back(i.propagated_distance);
    return vec;
}

double Secondaries::GetELost(const Geometry& geometry) const
{
    auto entry_point = GetEntryPoint(geometry);
    auto exit_point = GetExitPoint(geometry);
    return entry_point->energy - exit_point->energy;
}

std::unique_ptr<DynamicData> Secondaries::GetEntryPoint(
        const Geometry& geometry) const
{
    auto pos_0 = secondaries_.front().position;
    auto dir_0 = secondaries_.front().direction;
    if (geometry.IsEntering(pos_0, dir_0))
        return std::make_unique<DynamicData>(secondaries_.front());
    if (geometry.IsInside(pos_0, dir_0))
        return nullptr; // track starts in geometry

    for (unsigned int i = 0; i < secondaries_.size() - 1; i++) {
        auto pos_i = secondaries_.at(i).position;
        auto pos_f = secondaries_.at(i+1).position;
        auto dir_i = secondaries_.at(i).direction;
        auto dist_i_f = (pos_f - pos_i).magnitude();

        auto distance = geometry.DistanceToBorder(pos_i, dir_i).first;
        if (distance <= dist_i_f and distance >= 0) {
            if ( dist_i_f - distance < GEOMETRY_PRECISION)
                return std::make_unique<DynamicData>(secondaries_.at(i+1));
            auto entry_point = RePropagate(secondaries_.at(i), dir_i, distance);
            return std::make_unique<DynamicData>(entry_point);
        }
    }
    return nullptr; // No entry point found
}

std::unique_ptr<DynamicData> Secondaries::GetExitPoint(
        const Geometry &geometry) const
{
    auto pos_end = secondaries_.back().position;
    auto dir_end = secondaries_.back().direction;
    if (geometry.IsLeaving(pos_end, dir_end))
        return std::make_unique<DynamicData>(secondaries_.back());
    if (geometry.IsInside(pos_end, dir_end))
        return nullptr; // track ends inside geometry

    for (auto i = secondaries_.size() - 1; i > 0; i--) {
        auto pos_i = secondaries_.at(i-1).position;
        auto pos_f = secondaries_.at(i).position;
        auto dir_i = secondaries_.at(i-1).direction;
        auto dist_i_f = (pos_f - pos_i).magnitude();

        auto distance = geometry.DistanceToBorder(pos_f, -dir_i).first;
        if (distance <= dist_i_f and distance >= 0) {
            if ( dist_i_f - distance < GEOMETRY_PRECISION)
                return std::make_unique<DynamicData>(secondaries_.at(i-1));
            auto exit_point = RePropagate(secondaries_.at(i-1), dir_i,
                                          dist_i_f - distance);
            return std::make_unique<DynamicData>(exit_point);
        }
    }

    auto pos_0 = secondaries_.front().position;
    auto dir_0 = secondaries_.front().direction;
    if (geometry.IsLeaving(pos_0, dir_0))
        return std::make_unique<DynamicData>(secondaries_.front());

    return nullptr; // No exit point found
}

std::unique_ptr<DynamicData> Secondaries::GetClosestApproachPoint(const Geometry& geometry) const
{
    for (unsigned int i = 0; i < secondaries_.size(); i++) {
        auto sec_pos = secondaries_.at(i).position;
        auto sec_dir = secondaries_.at(i).direction;
        if (geometry.DistanceToClosestApproach(sec_pos, sec_dir) <= 0.) {
            if(std::abs(geometry.DistanceToClosestApproach(sec_pos, sec_dir))
                        < PARTICLE_POSITION_RESOLUTION)
                return std::make_unique<DynamicData>(secondaries_.at(i));
            if (i == 0)
                return std::make_unique<DynamicData>(secondaries_.front());
            auto prev_pos = secondaries_.at(i-1).position;
            auto prev_dir = secondaries_.at(i-1).direction;
            auto direction = sec_pos - prev_pos;
            direction.normalise();
            auto displacement = geometry.DistanceToClosestApproach(prev_pos,
                                                                   prev_dir);
            auto closest_approach = RePropagate(secondaries_.at(i-1), direction,
                                                displacement);
            return std::make_unique<DynamicData>(closest_approach);
        }
    }
    return std::make_unique<DynamicData>(secondaries_.back());
}

DynamicData Secondaries::RePropagate(const DynamicData &init,
                                     const Vector3D& direction,
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
    auto advance_grammage = density->Calculate(init.position, direction,
                                               displacement);
    auto E_f = utility.EnergyDistance(init.energy, advance_grammage);
    auto new_time = init.time + utility.TimeElapsed(
            init.energy, E_f,displacement, density->Evaluate(init.position));
    auto new_position = init.position + direction * displacement;
    auto new_propagated_distance = init.propagated_distance + displacement;
    return DynamicData((ParticleType)primary_def_->particle_type, new_position,
                       direction, E_f, new_time, new_propagated_distance);
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
