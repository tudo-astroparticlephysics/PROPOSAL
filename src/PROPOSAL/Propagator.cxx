#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/Secondaries.h"
#include "PROPOSAL/crosssection/ParticleDefaultCrossSectionList.h"
#include "PROPOSAL/crosssection/parametrization/Annihilation.h"
#include "PROPOSAL/crosssection/parametrization/Bremsstrahlung.h"
#include "PROPOSAL/crosssection/parametrization/Compton.h"
#include "PROPOSAL/crosssection/parametrization/EpairProduction.h"
#include "PROPOSAL/crosssection/parametrization/Ionization.h"
#include "PROPOSAL/crosssection/parametrization/MupairProduction.h"
#include "PROPOSAL/crosssection/parametrization/PhotoPairProduction.h"
#include "PROPOSAL/crosssection/parametrization/PhotoQ2Integration.h"
#include "PROPOSAL/crosssection/parametrization/PhotoRealPhotonAssumption.h"
#include "PROPOSAL/crosssection/parametrization/WeakInteraction.h"
#include "PROPOSAL/geometry/GeometryFactory.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/propagation_utility/InteractionBuilder.h"
#include "PROPOSAL/propagation_utility/TimeBuilder.h"
#include "PROPOSAL/propagation_utility/DecayBuilder.h"
#include "PROPOSAL/propagation_utility/ContRandBuilder.h"
#include "PROPOSAL/scattering/ScatteringFactory.h"
#include "PROPOSAL/density_distr/density_distr.h"
#include <fstream>
#include "PROPOSAL/Logging.h"

#include <iomanip>

using namespace PROPOSAL;
using std::make_shared;
using std::string;

Propagator::Propagator(const ParticleDef& p_def, std::vector<Sector> sectors)
    : p_def(p_def)
    , sector_list(sectors)
{
}

Secondaries Propagator::Propagate(
    const DynamicData& initial_particle, double max_distance, double min_energy)
{
    Secondaries track(std::make_shared<ParticleDef>(p_def), sector_list);
    track.push_back(initial_particle);
    track.back().SetType(p_def.particle_type);

    auto current_sector = GetCurrentSector(
            track.back().GetPosition(), track.back().GetDirection());
    auto rnd
        = std::bind(&RandomGenerator::RandomDouble, &RandomGenerator::Get());

    int advancement_type;
    auto continue_propagation = true;

    // TODO: How to get accurate low information?
    auto InteractionEnergy
        = std::array<double, 3>{std::max(min_energy, p_def.mass), 0., 0.};
    while (continue_propagation) {
        //std::cout << "E: " << track.back().GetEnergy() << ", d: "  << track.back().GetPropagatedDistance() << std::endl;
        auto &utility = get<UTILITY>(current_sector);
        auto &density = get<DENSITY_DISTR>(current_sector);

        InteractionEnergy[Decay] = utility.EnergyDecay(track.back().GetEnergy(),
                                                       rnd, density->Evaluate(track.back().GetPosition()));
        InteractionEnergy[Stochastic]
                = utility.EnergyInteraction(track.back().GetEnergy(), rnd);

        //std::cout << "Decay: " << InteractionEnergy[Decay] << ", " << "Stochastic: " << InteractionEnergy[Stochastic] << std::endl;

        auto next_interaction_type = maximize(InteractionEnergy);
        auto energy_at_next_interaction = InteractionEnergy[next_interaction_type];

        track.push_back(track.back());
        advancement_type = AdvanceParticle(track.back(),
                                           energy_at_next_interaction,
                                           max_distance, rnd, current_sector);
        switch (advancement_type) {
            case ReachedInteraction :
                switch (next_interaction_type) {
                    case Stochastic:
                        track.push_back(track.back());
                        DoStochasticInteraction(track.back(), utility, rnd);
                        if (track.back().GetEnergy() <= InteractionEnergy[MinimalE])
                            continue_propagation = false;
                        break;
                    case Decay:
                        track.push_back(track.back());
                        track.back().SetType(InteractionType::Decay);
                        continue_propagation = false;
                        break;
                    case MinimalE:
                        continue_propagation = false;
                        break;
                }
                break;
            case ReachedBorder :
                current_sector = GetCurrentSector(track.back().GetPosition(), track.back().GetDirection());
                break;
            case ReachedMaxDistance :
                continue_propagation = false;
                break;
        }
    }
    return track;
}

void Propagator::DoStochasticInteraction(DynamicData& p_cond,
    PropagationUtility& utility, std::function<double()> rnd)
{
    InteractionType loss_type;
    std::shared_ptr<const Component> comp;
    double loss_energy;
    std::tie(loss_type, comp, loss_energy)
        = utility.EnergyStochasticloss(p_cond.GetEnergy(), rnd());

    /* auto deflection_angles = stochastic_cross->StochasticDeflection( */
    /*     p_cond.GetEnergy(), loss_energy); // TODO: ugly */
    /* p_cond.DeflectDirection( */
    /*     get<0>(deflection_angles), get<1>(deflection_angles)); */

    p_cond.SetEnergy(p_cond.GetEnergy() - loss_energy);
    p_cond.SetType(loss_type);
}

int Propagator::AdvanceParticle(DynamicData& p_cond, double E_f,
        double max_distance, std::function<double()> rnd, Sector& current_sector)
{
    assert(max_distance > 0);
    assert(E_f >= 0);

    // TODO: For NoScattering, these random numbers are not used
    auto rnd_scattering = std::array<double, 4>{ rnd(), rnd(), rnd(), rnd() };

    auto& utility = get<UTILITY>(current_sector);
    auto& density = get<DENSITY_DISTR>(current_sector);
    auto& geometry = get<GEOMETRY>(current_sector);

    auto grammage_next_interaction = utility.LengthContinuous(p_cond.GetEnergy(), E_f);
    auto max_distance_left = max_distance - p_cond.GetPropagatedDistance();
    assert(max_distance_left > 0);

    Vector3D mean_direction, new_direction;
    std::tie(mean_direction, new_direction) = utility.DirectionsScatter(
            grammage_next_interaction, p_cond.GetEnergy(), E_f, p_cond.GetDirection(),
        rnd_scattering);

    double distance_next_interaction;
    try {
        distance_next_interaction = density->Correct(p_cond.GetPosition(),
                  mean_direction, grammage_next_interaction, max_distance_left);
    } catch (const DensityException&) {
        Logging::Get("proposal.propagator")->debug("Interaction point exceeds "
                     "maximum propagation distance or lies in infinity.");
        distance_next_interaction = INF;
    }
    auto new_position = p_cond.GetPosition() + distance_next_interaction * mean_direction;

    auto AdvanceDistance = std::array<double, 3>{0., 0., 0.};
    AdvanceDistance[ReachedInteraction] = distance_next_interaction;
    AdvanceDistance[ReachedMaxDistance] = max_distance_left;
    AdvanceDistance[ReachedBorder] = CalculateDistanceToBorder(p_cond.GetPosition(), mean_direction, *geometry);

    int advancement_type = minimize(AdvanceDistance);
    double advance_distance = AdvanceDistance[advancement_type];

    if(advancement_type != ReachedInteraction) {
        double advance_grammage;
        double control_distance;
        //std::cout << "AdvanceParticle can't reach interaction, varying propagation step..." << std::endl;
        do {
            advance_distance = AdvanceDistance[advancement_type];
            advance_grammage = density->Calculate(p_cond.GetPosition(),
                                       p_cond.GetDirection(), advance_distance);
            E_f = utility.EnergyDistance(p_cond.GetEnergy(), advance_grammage);

            std::tie(mean_direction, new_direction) = utility.DirectionsScatter(
                    advance_grammage, p_cond.GetEnergy(), E_f,
                    p_cond.GetDirection(), rnd_scattering);

            try {
                AdvanceDistance[ReachedInteraction] = density->Correct(
                        p_cond.GetPosition(), mean_direction,
                        grammage_next_interaction, max_distance_left);
            } catch (const DensityException&) {
                AdvanceDistance[ReachedInteraction] = INF;
            }

            AdvanceDistance[ReachedBorder] = CalculateDistanceToBorder(
                    p_cond.GetPosition(), mean_direction, *geometry);
            advancement_type = minimize(AdvanceDistance);
            control_distance = AdvanceDistance[advancement_type];
            //std::cout << "Step - old_distance: " << advance_distance << ", new distance: " << control_distance << ", advancement type " << advancement_type << std::endl;
        } while (std::abs(advance_distance - control_distance)
                 > PARTICLE_POSITION_RESOLUTION);
        //std::cout << "Difference negligible, use control_distance" << std::endl;
        advance_distance = control_distance;
        new_position = p_cond.GetPosition() + advance_distance * mean_direction;
    }

    p_cond.SetTime(p_cond.GetTime() + utility.TimeElapsed(p_cond.GetEnergy(),
            E_f, advance_distance, density->Evaluate(p_cond.GetPosition())));
    p_cond.SetPosition(new_position);
    p_cond.SetDirection(new_direction);
    p_cond.SetPropagatedDistance(p_cond.GetPropagatedDistance()
        + advance_distance);
    p_cond.SetEnergy(utility.EnergyRandomize(p_cond.GetEnergy(), E_f, rnd));
    p_cond.SetType(InteractionType::ContinuousEnergyLoss);

    return advancement_type;
}

double Propagator::CalculateDistanceToBorder(const Vector3D& position,
    const Vector3D& direction, const Geometry& current_geometry)
{
    auto distance_border
        = current_geometry.DistanceToBorder(position, direction).first;
    double tmp_distance;
    for (auto& sector : sector_list) {
        auto& geometry = get<GEOMETRY>(sector);
        if (geometry->GetHierarchy() > current_geometry.GetHierarchy()) {
            tmp_distance = geometry->DistanceToBorder(position, direction).first;
            if(tmp_distance >= 0)
                distance_border = std::min(distance_border, tmp_distance);
        }
    }
    return distance_border;
}

int Propagator::maximize(const std::array<double, 3>& InteractionEnergies)
{
    auto max_element_ref = std::max_element(
        InteractionEnergies.begin(), InteractionEnergies.end());
    return std::distance(InteractionEnergies.begin(), max_element_ref);
}

int Propagator::minimize(const std::array<double, 3>& AdvanceDistances)
{
    auto min_element_ref = std::min_element(
            AdvanceDistances.begin(), AdvanceDistances.end());
    return std::distance(AdvanceDistances.begin(), min_element_ref);
}


Sector Propagator::GetCurrentSector(
    const Vector3D& position, const Vector3D& direction)
{
    auto potential_sec = std::vector<Sector*>{};
    for (auto& sector : sector_list) {
        if (get<GEOMETRY>(sector)->IsInside(position, direction))
            potential_sec.push_back(&sector);
    }

    if (potential_sec.empty())
        throw std::logic_error(
            "Propagator: No sector defined at current particle position.");

    auto highest_sector_iter = std::max_element(
        potential_sec.begin(), potential_sec.end(), [](Sector* a, Sector* b) {
            return get<GEOMETRY>(*a)->GetHierarchy()
                < get<GEOMETRY>(*b)->GetHierarchy();
        });

    return **highest_sector_iter;
}

// Init methods

nlohmann::json Propagator::ParseConfig(const string& config_file)
{
    nlohmann::json json_config;
    try {
        string expanded_config_file_path
            = Helper::ResolvePath(config_file, true);
        std::ifstream infilestream(expanded_config_file_path);
        infilestream >> json_config;
    } catch (const nlohmann::json::parse_error& e) {
        Logging::Get("proposal.propagator")->critical("Unable parse \"%s\" as json file", config_file.c_str());
    }
    return json_config;
}

Propagator::GlobalSettings::GlobalSettings(const nlohmann::json& config_global)
{
    if (config_global.contains("medium"))
        medium = CreateMedium(config_global["medium"].get<std::string>());
    cross = config_global.value("CrossSections", "");
    if (config_global.contains("cuts"))
        cuts = make_shared<EnergyCutSettings>(config_global["cuts"]);
    do_exact_time = config_global.value("exact_time", true);
    do_interpolation = config_global.value("do_interpolation", true);
    scattering = config_global.value("scattering", "NoScattering");
}

Propagator::GlobalSettings::GlobalSettings()
{
    medium = nullptr;
    cross = nullptr;
    do_exact_time = true;
    do_interpolation = true;
    scattering = "NoScattering";
}
