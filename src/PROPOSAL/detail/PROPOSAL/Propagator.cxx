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
#include "PROPOSAL/crosssection/Factories/PhotoeffectFactory.h"
#include "PROPOSAL/crosssection/Factories/PhotoMuPairProductionFactory.h"
#include "PROPOSAL/crosssection/Factories/PhotoPairProductionFactory.h"
#include "PROPOSAL/crosssection/Factories/PhotonuclearFactory.h"
#include "PROPOSAL/crosssection/Factories/PhotoproductionFactory.h"
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

        advancement_type = AdvanceParticle(
                state, energy_at_next_interaction, max_distance, rnd,
                current_sector, next_interaction_type == MinimalE,
                InteractionEnergy[MinimalE]);

        // If the particle is on the sector border before the continuous step is
        // performed in 'AdvanceParticle', we might enter a different sector due
        // to multiple scattering. Therefore, we have to update current_sector.
        utility = get<UTILITY>(current_sector);
        density = get<DENSITY_DISTR>(current_sector);

        track.push_back(state, InteractionType::ContinuousEnergyLoss);

        switch (advancement_type) {
        case ReachedInteraction:
            switch (next_interaction_type) {
            case Stochastic: {
                auto loss = DoStochasticInteraction(state, utility, rnd);
                if (loss.type != InteractionType::Undefined)
                    track.push_back(state, loss.type, loss.comp_hash);
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

Interaction::Loss Propagator::DoStochasticInteraction(ParticleState& p_cond,
    PropagationUtility& utility, std::function<double()> rnd)
{
    auto loss = utility.EnergyStochasticloss(p_cond.energy, rnd());

    p_cond.direction = utility.DirectionDeflect(loss.type, p_cond.energy,
        p_cond.energy * (1. - loss.v_loss), p_cond.direction, rnd, loss.comp_hash);
    p_cond.energy = p_cond.energy * (1. - loss.v_loss);

    return loss;
}

int Propagator::AdvanceParticle(ParticleState &state,
    const double energy_next_interaction, const double final_distance,
    std::function<double()> rnd_generator, Sector& current_sector,
    bool min_energy_step, const double min_energy) {

    auto& utility = get<UTILITY>(current_sector);
    auto& density = get<DENSITY_DISTR>(current_sector);
    auto& geometry = get<GEOMETRY>(current_sector);

    double energy = energy_next_interaction; // final energy of proposed step
    double grammage = -1; // grammage of proposed step
    double distance = -1; // geometrical distance of proposed step

    // Calculate maximal allowed length of step (limit due to final_distance)
    const double max_distance = final_distance - state.propagated_distance;

    // Calculate grammage until next stochastic interaction
    double grammage_next_interaction = utility.LengthContinuous(
            state.energy, energy_next_interaction);

    int advancement_type;
    Cartesian3D mean_direction, new_direction; // proposed scattering

    // This lambda expression ensures we always use the same 4 random numbers
    std::array<double, 4> random_numbers = {-1, -1, -1, -1};
    int i = 0;
    auto rnd = [&rnd_generator, &i, &random_numbers]() {
        if (random_numbers[i%4] == -1)
            random_numbers[i%4] = rnd_generator();
        return random_numbers[i++%4];
    };

    int num_steps = 0; // count number of iteration steps
    bool backscatter = false;

    // Iterate combinations of step lengths and scattering angles until we have
    // reached an interaction, a sector border or the maximal propagation distance
    do {
        num_steps++;
        // Calculate grammage, energy and distance for step
        if (energy != -1 && distance == -1) {
            // Calculate grammage and distance from given energy
            grammage = utility.LengthContinuous(state.energy, energy);
            try {
                distance = density->Correct(state.position, state.direction, grammage, max_distance);
            } catch (const DensityException&) {
                distance = INF;
            }
        } else if (energy == -1 && distance != -1) {
            // Calculate energy and grammage from given distance
            auto grammage_step = density->Calculate(state.position, state.direction, distance);
            if (grammage_step < grammage_next_interaction) {
                grammage = grammage_step;
                energy = utility.EnergyDistance(state.energy, grammage);
            } else {
                // we are unable to reach `distance` before we reach the next interaction
                // this means we are stuck in a loop, and need to discard the current set of random numbers
                for (auto& r: random_numbers) {
                    r = rnd_generator();
                }
                Logging::Get("proposal.propagator")->debug("Unable to find a valid combination of propagation step "
                                                           "length and multiple scattering angle for this set of "
                                                           "random numbers. Resample set of random numbers.");
                grammage = grammage_next_interaction;
                energy = energy_next_interaction;
                try {
                    distance = density->Correct(state.position, state.direction, grammage, max_distance);
                } catch (const DensityException&) {
                    distance = INF;
                }
            }
        } else {
            throw std::logic_error("Error in AdvanceParticle: Either both distance and final energy for the next "
                                   "iteration step are known, or neither of them are known. This should never happen "
                                   "and would indicate an algorithmic error!");
        }

        // Calculate scattering proposal
        std::tie(mean_direction, new_direction) = utility.DirectionsScatter(
                grammage, state.energy, energy, state.direction, rnd);

        // Check step
        double distance_to_border = CalculateDistanceToBorder(state.position, mean_direction, *geometry);
        bool is_inside = geometry->IsInside(state.position, mean_direction);

        if (num_steps > PropagationSettings::ADVANCE_PARTICLE_MAX_STEPS) {
            // too many iteration steps!
            Logging::Get("proposal.propagator")->warn("Maximal number of iteration step exceeded ({}). "
                                                      "Proposed propagation step is {} cm, while distance to border is "
                                                      "{} cm (difference of {} cm). Initial energy particle {} MeV, "
                                                      "particle energy after step {} MeV. Propagate to border and "
                                                      "continue propagation.",
                                                      PropagationSettings::ADVANCE_PARTICLE_MAX_STEPS, distance,
                                                      distance_to_border, std::abs(distance - distance_to_border),
                                                      state.energy, energy);
            distance = distance_to_border;
            grammage = density->Calculate(state.position, state.direction, distance);
            energy = utility.EnergyDistance(state.energy, grammage);
            advancement_type = ReachedBorder;
        } else if (!is_inside) {
            // Special case: We are on the sector border, but scattering back outside the current sector!
            // Update sector and recalculate values
            advancement_type = InvalidStep;
            auto new_sector = GetCurrentSector(state.position, mean_direction);
            utility = get<UTILITY>(new_sector);
            density = get<DENSITY_DISTR>(new_sector);
            geometry = get<GEOMETRY>(new_sector);
            grammage_next_interaction = utility.LengthContinuous(state.energy, energy_next_interaction);
            energy = energy_next_interaction;
            distance = -1;
            grammage = -1;

            // if we get in this case, this might mean that we are stuck in a loop.
            // this can happen if we backscatter in both the old and the new sector (if their medium is different).
            // sample new random numbers to avoid this.
            if (backscatter == true) {
                for (auto& r: random_numbers) {
                    r = rnd_generator();
                }
            }
            backscatter = true;
        } else if (distance <= distance_to_border && distance <= max_distance && energy == energy_next_interaction) {
            // reached interaction
            advancement_type = ReachedInteraction;
        } else if (distance <= distance_to_border && distance == max_distance) {
            // reached max distance
            advancement_type = ReachedMaxDistance;
        } else if (std::abs(distance - distance_to_border) <= PARTICLE_POSITION_RESOLUTION) {
            // reached geometry border
            advancement_type = ReachedBorder;
            distance = distance_to_border;
        } else {
            // iteration not finished! discard energy and grammage
            advancement_type = InvalidStep;
            distance = std::min(distance_to_border, max_distance);
            energy = -1;
            grammage = -1;
        }
    } while (advancement_type == InvalidStep);

    state.time = state.time + utility.TimeElapsed(state.energy, energy, grammage, density->Evaluate(state.position)); // TODO: should the energy passed here be the randomized energy or not?
    state.position = state.position + distance * mean_direction;
    state.direction = new_direction;
    state.propagated_distance = state.propagated_distance + distance;
    if (min_energy_step && advancement_type == ReachedInteraction)
        state.energy = energy; // we reached a specific energy, no randomization
    else
        state.energy = utility.EnergyRandomize(state.energy, energy, rnd, min_energy);

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
    std::ifstream f(config_file.c_str());
    if (!f.good()) {
        throw std::invalid_argument("No configuration file found under "
                                    "path "+ config_file);
    }
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
        = make_interaction(def.displacement_calc, crosss, do_interpol, false);
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
    if (config.contains("photoeffect"))
        cross.emplace_back(make_photoeffect(
            p_def, medium, config["photoeffect"]));
    if (config.contains("photomupair"))
        cross.emplace_back(make_photomupairproduction(
            p_def, medium, interpolate, config["photomupair"]));
    if (config.contains("photoproduction"))
        cross.emplace_back(make_photoproduction(
            p_def, medium, config["photoproduction"]));
    if (config.contains("photopair"))
        cross.emplace_back(make_photopairproduction(
                p_def, medium, interpolate, config["photopair"],
                density_correction));
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
