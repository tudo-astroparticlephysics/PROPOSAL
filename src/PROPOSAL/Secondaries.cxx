
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

void Secondaries::emplace_back(const int& type, const Vector3D& position,
    const Vector3D& direction, const double& energy,
    const double& parent_particle_energy, const double& time,
    const double& distance)
{
    secondaries_.emplace_back(type, position, direction, energy,
        parent_particle_energy, time, distance);
}
void Secondaries::emplace_back(const int& type)
{
    secondaries_.emplace_back(type);
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
        if (interaction_type == i.GetType())
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
        if (geometry.IsInside(i.GetPosition(), i.GetDirection()))
            sec.push_back(i);
    }
    return sec;
}

void Secondaries::DoDecay()
{
    for (auto it = secondaries_.begin(); it != secondaries_.end();) {
        if (it->GetType() == static_cast<int>(InteractionType::Decay)) {
            DynamicData decaying_particle(primary_def_->particle_type,
                it->GetPosition(), it->GetDirection(), it->GetEnergy(),
                it->GetParentParticleEnergy(), it->GetTime(),
                it->GetPropagatedDistance());
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
        vec.emplace_back(i.GetPosition());
    return vec;
}

std::vector<Vector3D> Secondaries::GetDirection() const
{
    std::vector<Vector3D> vec;
    for (auto i : secondaries_)
        vec.emplace_back(i.GetDirection());
    return vec;
}

std::vector<double> Secondaries::GetEnergy() const
{
    std::vector<double> vec;
    for (auto i : secondaries_)
        vec.emplace_back(i.GetEnergy());
    return vec;
}

std::vector<double> Secondaries::GetParentParticleEnergy() const
{
    std::vector<double> vec;
    for (auto i : secondaries_)
        vec.emplace_back(i.GetParentParticleEnergy());
    return vec;
}

std::vector<double> Secondaries::GetTime() const
{
    std::vector<double> vec;
    for (auto i : secondaries_)
        vec.emplace_back(i.GetTime());
    return vec;
}

std::vector<double> Secondaries::GetPropagatedDistance() const
{
    std::vector<double> vec;
    for (auto i : secondaries_)
        vec.emplace_back(i.GetPropagatedDistance());
    return vec;
}

double Secondaries::GetELost(const Geometry& geometry) const
{
    auto entry_point = GetEntryPoint(geometry);
    auto exit_point = GetExitPoint(geometry);
    return entry_point->GetEnergy() - exit_point->GetEnergy();
}

std::unique_ptr<DynamicData> Secondaries::GetEntryPoint(
        const Geometry& geometry) const
{
    auto pos_0 = secondaries_.front().GetPosition();
    auto dir_0 = secondaries_.front().GetDirection();
    if (geometry.IsEntering(pos_0, dir_0))
        return std::make_unique<DynamicData>(secondaries_.front());
    if (geometry.IsInside(pos_0, dir_0))
        return nullptr; // track starts in geometry

    for (unsigned int i = 0; i < secondaries_.size() - 1; i++) {
        auto pos_i = secondaries_.at(i).GetPosition();
        auto pos_f = secondaries_.at(i+1).GetPosition();
        auto dir_i = secondaries_.at(i).GetDirection();
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
    auto pos_end = secondaries_.back().GetPosition();
    auto dir_end = secondaries_.back().GetDirection();
    if (geometry.IsLeaving(pos_end, dir_end))
        return std::make_unique<DynamicData>(secondaries_.back());
    if (geometry.IsInside(pos_end, dir_end))
        return nullptr; // track ends inside geometry

    for (auto i = secondaries_.size() - 1; i > 0; i--) {
        auto pos_i = secondaries_.at(i-1).GetPosition();
        auto pos_f = secondaries_.at(i).GetPosition();
        auto dir_i = secondaries_.at(i-1).GetDirection();
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

    auto pos_0 = secondaries_.front().GetPosition();
    auto dir_0 = secondaries_.front().GetDirection();
    if (geometry.IsLeaving(pos_0, dir_0))
        return std::make_unique<DynamicData>(secondaries_.front());

    return nullptr; // No exit point found
}

std::unique_ptr<DynamicData> Secondaries::GetClosestApproachPoint(const Geometry& geometry) const
{
    for (unsigned int i = 0; i < secondaries_.size(); i++) {
        auto sec_pos = secondaries_.at(i).GetPosition();
        auto sec_dir = secondaries_.at(i).GetDirection();
        if (geometry.DistanceToClosestApproach(sec_pos, sec_dir) <= 0.) {
            if(std::abs(geometry.DistanceToClosestApproach(sec_pos, sec_dir))
                        < PARTICLE_POSITION_RESOLUTION)
                return std::make_unique<DynamicData>(secondaries_.at(i));
            if (i == 0)
                return std::make_unique<DynamicData>(secondaries_.front());
            auto prev_pos = secondaries_.at(i-1).GetPosition();
            auto prev_dir = secondaries_.at(i-1).GetDirection();
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
    auto current_sector = GetCurrentSector(init.GetPosition(),
                                           init.GetDirection());
    auto& utility = get<Propagator::UTILITY>(current_sector);
    auto& density = get<Propagator::DENSITY_DISTR>(current_sector);

    if (utility.collection.cont_rand != nullptr) {
        Logging::Get("proposal.secondaries")->warn("Calcuating specific points "
                     "for sectors where continuous randomization is activated "
                     "can lead to inconsistent results!");
    }
    auto advance_grammage = density->Calculate(init.GetPosition(), direction,
                                               displacement);
    auto E_f = utility.EnergyDistance(init.GetEnergy(), advance_grammage);
    auto new_time = utility.TimeElapsed(init.GetEnergy(), E_f, displacement,
                                        density->Evaluate(init.GetPosition()));
    auto new_position = init.GetPosition() + direction * displacement;
    auto new_propagated_distance = init.GetPropagatedDistance() + displacement;
    return DynamicData((int)InteractionType::ContinuousEnergyLoss, new_position,
                       direction, E_f, init.GetParentParticleEnergy(), new_time,
                       new_propagated_distance);
}

Sector Secondaries::GetCurrentSector(const Vector3D& position,
                                     const Vector3D& direction) const
{
    //TODO: this is essentially a dublicate of Propagator::GetCurrentSector
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
