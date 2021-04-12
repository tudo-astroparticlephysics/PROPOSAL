#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/Secondaries.h"
#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/Factories/AnnihilationFactory.h"
#include "PROPOSAL/crosssection/Factories/BremsstrahlungFactory.h"
#include "PROPOSAL/crosssection/Factories/ComptonFactory.h"
#include "PROPOSAL/crosssection/Factories/EpairProductionFactory.h"
#include "PROPOSAL/crosssection/Factories/IonizationFactory.h"
#include "PROPOSAL/crosssection/Factories/MupairProductionFactory.h"
#include "PROPOSAL/crosssection/Factories/PhotoPairProductionFactory.h"
#include "PROPOSAL/crosssection/Factories/PhotonuclearFactory.h"
#include "PROPOSAL/crosssection/Factories/WeakInteractionFactory.h"
#include "PROPOSAL/crosssection/ParticleDefaultCrossSectionList.h"
#include "PROPOSAL/density_distr/density_distr.h"
#include "PROPOSAL/geometry/GeometryFactory.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/propagation_utility/ContRandBuilder.h"
#include "PROPOSAL/propagation_utility/DecayBuilder.h"
#include "PROPOSAL/propagation_utility/InteractionBuilder.h"
#include "PROPOSAL/propagation_utility/TimeBuilder.h"
#include "PROPOSAL/scattering/ScatteringFactory.h"
#include <fstream>

#include <iomanip>

using namespace PROPOSAL;
using std::get;
using std::string;

Propagator::Propagator(const ParticleDef& p_def, std::vector<Sector> sectors)
    : p_def(p_def)
    , sector_list(sectors)
{
}

Propagator::Propagator(const ParticleDef& p_def, const nlohmann::json& config)
    : p_def(p_def)
{
    GlobalSettings global;
    if (config.contains("global"))
        global = GlobalSettings(config["global"]);
    if (config.contains("sectors")) {
        assert(config["sectors"].is_array());
        for (const auto& json_sector : config.at("sectors")) {
            InitializeSectorFromJSON(p_def, json_sector, global);
        }
    } else {
        throw std::invalid_argument("No sector array found in json object");
    }
}

Secondaries Propagator::Propagate(const ParticleState& initial_particle,
    double max_distance, double min_energy, unsigned int hierarchy_condition)
{
    Secondaries track(std::make_shared<ParticleDef>(p_def), sector_list);

    track.push_back(initial_particle, InteractionType::ContinuousEnergyLoss);
    auto state = ParticleState(initial_particle);

    auto current_sector = GetCurrentSector(state.position, state.direction);
    auto rnd
        = std::bind(&RandomGenerator::RandomDouble, &RandomGenerator::Get());

    int advancement_type;
    auto continue_propagation = true;

    std::array<double, 3> InteractionEnergy;
    while (continue_propagation) {
        auto& utility = get<UTILITY>(current_sector);
        auto& density = get<DENSITY_DISTR>(current_sector);

        InteractionEnergy[MinimalE] = std::max(
                min_energy, utility.collection.displacement_calc->GetLowerLim());
        InteractionEnergy[Decay] = utility.EnergyDecay(
            state.energy, rnd, density->Evaluate(state.position));
        InteractionEnergy[Stochastic]
            = utility.EnergyInteraction(state.energy, rnd);

        auto next_interaction_type = maximize(InteractionEnergy);
        auto energy_at_next_interaction
            = InteractionEnergy[next_interaction_type];

        advancement_type = AdvanceParticle(state, energy_at_next_interaction,
            max_distance, rnd, current_sector);
        track.push_back(state, InteractionType::ContinuousEnergyLoss);

        switch (advancement_type) {
        case ReachedInteraction:
            switch (next_interaction_type) {
            case Stochastic: {
                auto type = DoStochasticInteraction(state, utility, rnd);
                track.push_back(state, type);
                if (state.energy <= InteractionEnergy[MinimalE])
                    continue_propagation = false;
                break;
            }
            case Decay: {
                track.push_back(state, InteractionType::Decay);
                continue_propagation = false;
                break;
            }
            case MinimalE: {
                continue_propagation = false;
                break;
            }
            }
            break;
        case ReachedBorder: {
            auto hierarchy_i = get<GEOMETRY>(current_sector)->GetHierarchy();
            current_sector = GetCurrentSector(state.position, state.direction);
            auto hierarchy_f = get<GEOMETRY>(current_sector)->GetHierarchy();
            if (hierarchy_i > hierarchy_condition
                && hierarchy_f < hierarchy_condition)
                continue_propagation = false;
            break;
        }
        case ReachedMaxDistance:
            continue_propagation = false;
            break;
        }
    }
    return track;
}

InteractionType Propagator::DoStochasticInteraction(ParticleState& p_cond,
    PropagationUtility& utility, std::function<double()> rnd)
{
    auto loss = utility.EnergyStochasticloss(p_cond.energy, rnd());

    p_cond.direction = utility.DirectionDeflect(loss.type, p_cond.energy,
        p_cond.energy * (1. - loss.v_loss), p_cond.direction, rnd);
    p_cond.energy = p_cond.energy * (1. - loss.v_loss);

    return loss.type;
}

int Propagator::AdvanceParticle(ParticleState& state, double E_f,
    double max_distance, std::function<double()> rnd, Sector& current_sector)
{
    assert(max_distance > 0);
    assert(E_f >= 0);

    auto& utility = get<UTILITY>(current_sector);
    auto& density = get<DENSITY_DISTR>(current_sector);
    auto& geometry = get<GEOMETRY>(current_sector);

    auto grammage_next_interaction
        = utility.LengthContinuous(state.energy, E_f);
    auto max_distance_left = max_distance - state.propagated_distance;
    assert(max_distance_left > 0);

    Cartesian3D mean_direction, new_direction;
    std::tie(mean_direction, new_direction) = utility.DirectionsScatter(
        grammage_next_interaction, state.energy, E_f, state.direction, rnd);

    double distance_next_interaction;
    try {
        distance_next_interaction = density->Correct(state.position,
            mean_direction, grammage_next_interaction, max_distance_left);
    } catch (const DensityException&) {
        Logging::Get("proposal.propagator")
            ->debug("Interaction point exceeds "
                    "maximum propagation distance or lies in infinity.");
        distance_next_interaction = INF;
    }
    auto new_position
        = state.position + distance_next_interaction * mean_direction;

    auto AdvanceDistance = std::array<double, 3> { 0., 0., 0. };
    AdvanceDistance[ReachedInteraction] = distance_next_interaction;
    AdvanceDistance[ReachedMaxDistance] = max_distance_left;
    AdvanceDistance[ReachedBorder]
        = CalculateDistanceToBorder(state.position, mean_direction, *geometry);

    int advancement_type = minimize(AdvanceDistance);
    double advance_distance = AdvanceDistance[advancement_type];
    double advance_grammage = grammage_next_interaction;

    if (advancement_type != ReachedInteraction) {
        double control_distance;
        do {
            advance_distance = AdvanceDistance[advancement_type];
            advance_grammage = density->Calculate(
                state.position, state.direction, advance_distance);
            E_f = utility.EnergyDistance(state.energy, advance_grammage);

            std::tie(mean_direction, new_direction) = utility.DirectionsScatter(
                advance_grammage, state.energy, E_f, state.direction, rnd);

            try {
                AdvanceDistance[ReachedInteraction]
                    = density->Correct(state.position, mean_direction,
                        grammage_next_interaction, max_distance_left);
            } catch (const DensityException&) {
                AdvanceDistance[ReachedInteraction] = INF;
            }

            AdvanceDistance[ReachedBorder] = CalculateDistanceToBorder(
                state.position, mean_direction, *geometry);
            advancement_type = minimize(AdvanceDistance);
            control_distance = AdvanceDistance[advancement_type];
        } while (std::abs(advance_distance - control_distance)
            > PARTICLE_POSITION_RESOLUTION);
        advance_distance = control_distance;
        advance_grammage = density->Calculate(
            state.position, mean_direction, advance_distance);
        new_position = state.position + advance_distance * mean_direction;
    }

    state.time = state.time
        + utility.TimeElapsed(state.energy, E_f, advance_grammage,
            density->Evaluate(state.position));
    state.position = new_position;
    state.direction = new_direction;
    state.propagated_distance = state.propagated_distance + advance_distance;
    state.energy = utility.EnergyRandomize(state.energy, E_f, rnd);

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
            tmp_distance
                = geometry->DistanceToBorder(position, direction).first;
            if (tmp_distance >= 0)
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
    auto min_element_ref
        = std::min_element(AdvanceDistances.begin(), AdvanceDistances.end());
    return std::distance(AdvanceDistances.begin(), min_element_ref);
}

Sector Propagator::GetCurrentSector(
    const Vector3D& position, const Vector3D& direction)
{
    auto potential_sec = std::vector<Sector*> {};
    for (auto& sector : sector_list) {
        if (get<GEOMETRY>(sector)->IsInside(position, direction))
            potential_sec.push_back(&sector);
    }

    if (potential_sec.empty()) {
        auto spherical_position = Cartesian3D(position);
        Logging::Get("proposal.propagator")->critical("No sector defined at particle position {}, {}, {}.",
                                                      spherical_position.GetX(),
                                                      spherical_position.GetY(),
                                                      spherical_position.GetZ());
    }
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
        if (config_file== "")
            throw std::invalid_argument("Resolved path empty.");
        std::ifstream infilestream(config_file);
        infilestream >> json_config;
    } catch (const nlohmann::json::parse_error& e) {
        Logging::Get("proposal.propagator")->critical("Unable to parse {} as "
            "json file", config_file.c_str());
        throw;
    } catch (const std::invalid_argument& e) {
        Logging::Get("proposal.propagator")->critical("Could not find json"
            " file in path {}.", config_file.c_str());
        throw;
    }
    return json_config;
}

void Propagator::InitializeSectorFromJSON(const ParticleDef& p_def,
    const nlohmann::json& json_sector, GlobalSettings global)
{
    bool do_interpolation
        = json_sector.value("do_interpolation", global.do_interpolation);
    bool do_exact_time = json_sector.value("exact_time", global.do_exact_time);
    auto scattering_config = json_sector.value("scattering", global.scattering);
    std::shared_ptr<Medium> medium = global.medium;
    if (json_sector.contains("medium")) {
        medium = CreateMedium(json_sector["medium"].get<std::string>());
    } else if (medium == nullptr) {
        throw std::invalid_argument(
            "Neither a specific Sector medium nor a global medium is defined.");
    }
    std::shared_ptr<EnergyCutSettings> cuts = global.cuts;
    if (json_sector.contains("cuts")) {
        cuts = std::make_shared<EnergyCutSettings>(
            EnergyCutSettings(json_sector["cuts"]));
    } else if (cuts == nullptr) {
        throw std::invalid_argument("Neither a specific Sector EnergyCut nor a "
                                    "global EnergyCut is defined.");
    }
    nlohmann::json density_distr = { { "type", "homogeneous" },
        { "mass_density", medium->GetMassDensity() } };
    if (json_sector.contains("density_distribution"))
        density_distr = json_sector["density_distribution"];

    auto cross_config = json_sector.value("CrossSections", global.cross);
    PropagationUtility::Collection collection;
    if (!cross_config.empty()) {
        double density_correction
            = density_distr.value("mass_density", medium->GetMassDensity());
        density_correction /= medium->GetMassDensity();
        auto crosss = CreateCrossSectionList(p_def, *medium, cuts,
            do_interpolation, density_correction, cross_config);
        collection = CreateUtility(crosss, medium, cuts->GetContRand(),
            do_interpolation, do_exact_time, scattering_config);
    } else {
        auto std_crosss
            = GetStdCrossSections(p_def, *medium, cuts, do_interpolation);
        collection = CreateUtility(std_crosss, medium, cuts->GetContRand(),
            do_interpolation, do_exact_time, scattering_config);
    }
    auto utility = PropagationUtility(collection);

    if (json_sector.contains("geometries")) {
        assert(json_sector["geometries"].is_array());
        for (const auto& json_geometry : json_sector.at("geometries")) {
            auto geometry = CreateGeometry(json_geometry);
            auto density = CreateDensityDistribution(density_distr);
            sector_list.emplace_back(
                std::make_tuple(geometry, utility, density));
        }
    } else {
        throw std::invalid_argument(
            "At least one geometry must be defined for each sector");
    }
}

PropagationUtility::Collection Propagator::CreateUtility(
    std::vector<std::shared_ptr<CrossSectionBase>> crosss,
    std::shared_ptr<Medium> medium, bool do_cont_rand, bool do_interpol,
    bool do_exact_time, nlohmann::json scatter)
{
    PropagationUtility::Collection def;
    def.displacement_calc = make_displacement(crosss, do_interpol);
    def.interaction_calc
        = make_interaction(def.displacement_calc, crosss, do_interpol);
    if (!scatter.empty())
        def.scattering
            = make_scattering(scatter, p_def, *medium, crosss, do_interpol);
    if (std::isfinite(p_def.lifetime))
        def.decay_calc = make_decay(crosss, p_def, do_interpol);
    if (do_cont_rand)
        def.cont_rand = make_contrand(crosss, do_interpol);
    if (do_exact_time) {
        def.time_calc = make_time(crosss, p_def, do_interpol);
    } else {
        def.time_calc = std::make_shared<ApproximateTimeBuilder>();
    }
    return def;
}

std::vector<std::shared_ptr<CrossSectionBase>>
Propagator::CreateCrossSectionList(const ParticleDef& p_def,
    const Medium& medium, std::shared_ptr<const EnergyCutSettings> cuts,
    bool interpolate, double density_correction, const nlohmann::json& config)
{
    std::vector<std::shared_ptr<CrossSectionBase>> cross;

    if (config.contains("annihilation"))
        cross.emplace_back(make_annihilation(
            p_def, medium, interpolate, config["annihilation"]));
    if (config.contains("brems"))
        cross.emplace_back(make_bremsstrahlung(p_def, medium, cuts, interpolate,
            config["brems"], density_correction));
    if (config.contains("compton"))
        cross.emplace_back(
            make_compton(p_def, medium, cuts, interpolate, config["compton"]));
    if (config.contains("epair"))
        cross.emplace_back(make_epairproduction(p_def, medium, cuts,
            interpolate, config["epair"], density_correction));
    if (config.contains("ioniz"))
        cross.emplace_back(
            make_ionization(p_def, medium, cuts, interpolate, config["ioniz"]));
    if (config.contains("mupair"))
        cross.emplace_back(make_mupairproduction(
            p_def, medium, cuts, interpolate, config["mupair"]));
    if (config.contains("photo")) {
        try {
            cross.emplace_back(make_photonuclearreal(
                p_def, medium, cuts, interpolate, config["photo"]));
        } catch (std::invalid_argument& e) {
            cross.emplace_back(make_photonuclearQ2(
                p_def, medium, cuts, interpolate, config["photo"]));
        }
    }
    if (config.contains("photopair"))
        cross.emplace_back(make_photopairproduction(
            p_def, medium, interpolate, config["photopair"]));
    if (config.contains("weak"))
        cross.emplace_back(
            make_weakinteraction(p_def, medium, interpolate, config["weak"]));
    return cross;
}

Propagator::GlobalSettings::GlobalSettings(const nlohmann::json& config_global)
{
    if (config_global.contains("medium"))
        medium = CreateMedium(config_global["medium"].get<std::string>());
    if (config_global.contains("CrossSections"))
        cross = config_global["CrossSections"];
    if (config_global.contains("cuts"))
        cuts = std::make_shared<EnergyCutSettings>(config_global["cuts"]);
    if (config_global.contains("scattering"))
        scattering = config_global["scattering"];
    do_exact_time = config_global.value("exact_time", true);
    do_interpolation = config_global.value("do_interpolation", true);
}

Propagator::GlobalSettings::GlobalSettings()
{
    medium = nullptr;
    cross = {};
    do_exact_time = true;
    do_interpolation = true;
    scattering = {};
}
