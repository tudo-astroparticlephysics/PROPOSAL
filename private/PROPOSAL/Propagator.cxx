#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/Secondaries.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/geometry/GeometryFactory.h"
#include "PROPOSAL/scattering/ScatteringFactory.h"
#include "PROPOSAL/crossection/ParticleDefaultCrossSectionList.h"
#include "PROPOSAL/crossection/factories/AnnihilationFactory.h"
#include "PROPOSAL/crossection/factories/BremsstrahlungFactory.h"
#include "PROPOSAL/crossection/factories/ComptonFactory.h"
#include "PROPOSAL/crossection/factories/EpairProductionFactory.h"
#include "PROPOSAL/crossection/factories/IonizationFactory.h"
#include "PROPOSAL/crossection/factories/MupairProductionFactory.h"
#include "PROPOSAL/crossection/factories/PhotonuclearFactory.h"
#include "PROPOSAL/crossection/factories/PhotoPairFactory.h"
#include "PROPOSAL/crossection/factories/WeakInteractionFactory.h"
#include <fstream>

using namespace PROPOSAL;

Propagator::Propagator(const ParticleDef& p_def, const std::string& config_file)
                       : Propagator(p_def, ParseConfig(config_file)){
}

Propagator::Propagator(const ParticleDef& p_def, const nlohmann::json& config)
                        : p_def(p_def){

    GlobalSettings global;

    if (config.contains("global")) {
        global = GlobalSettings(config["global"]);
    }

    if (config.contains("interpolation")) {
        interpol_def_global = std::make_shared<InterpolationDef>(config["interpolation"]);
    } else {
        interpol_def_global = std::make_shared<InterpolationDef>();
    }


    if (config.contains("sectors")) {
        assert(config["sectors"].is_array());
        for (const auto& json_sector : config.at("sectors")) {
            InitializeSectorFromJSON(p_def, json_sector, global);
        }
    } else {
        throw std::invalid_argument("No sector array found in json object");
    }
}

Propagator::Propagator(const ParticleDef& p_def, std::vector<Sector> sectors) : p_def(p_def), sector_list(sectors) {
}

std::vector<DynamicData> Propagator::Propagate(const DynamicData& initial_particle, double max_distance, double min_energy) {
    std::vector<DynamicData> track;
    track.push_back(initial_particle);
    track.back().SetType(p_def.particle_type);

    auto current_sector = ChooseCurrentSector(track.back().GetPosition(), track.back().GetDirection());
    std::function<double()> rnd = std::bind(&RandomGenerator::RandomDouble, &RandomGenerator::Get());

    int next_interaction_type;
    double energy_at_next_interaction;
    double distance_to_next_interaction;
    bool sector_changed;
    bool continue_propagation = true;
    min_energy = std::max(min_energy, p_def.mass); //TODO: How to get accurate low information?
    while(continue_propagation) {
        std::array<double, 3> InteractionEnergy;
        InteractionEnergy[Decay] = get<UTILITY>(current_sector).EnergyDecay(track.back().GetEnergy(), rnd); // DEBUG: rnd needs to be corrected by density (as a first approximation)
        InteractionEnergy[Stochastic] = get<UTILITY>(current_sector).EnergyInteraction(track.back().GetEnergy(), rnd);
        InteractionEnergy[MinimalE] = min_energy;

        next_interaction_type = maximize(InteractionEnergy);
        energy_at_next_interaction = InteractionEnergy[next_interaction_type];
        distance_to_next_interaction = get<UTILITY>(current_sector).LengthContinuous(track.back().GetEnergy(), energy_at_next_interaction, 0.); // DEBUG: return value is grammage

        double max_distance_left = max_distance - track.back().GetPropagatedDistance(); // DEBUG: we need to make grammage out of max_distance_left to compare it in the next step (using density_distr->integrate)
        if (max_distance_left < distance_to_next_interaction) {
            distance_to_next_interaction = max_distance_left; // DEBUG: This is grammage now
            next_interaction_type = MaxDistance;
            energy_at_next_interaction = get<UTILITY>(current_sector).EnergyDistance(track.back().GetEnergy(), distance_to_next_interaction); // DEBUG: Insert grammage here
        }

        double approaching_sector_distance = std::numeric_limits<double>::infinity(); //TODO: Add new concept of smaller steps if approaching sector
        if (approaching_sector_distance < distance_to_next_interaction) { // DEBUG: both needs to be grammage
            distance_to_next_interaction = approaching_sector_distance; // DEBUG: This is grammage now
            next_interaction_type = ApproachingSector;
            energy_at_next_interaction = get<UTILITY>(current_sector).EnergyDistance(track.back().GetEnergy(), distance_to_next_interaction); // DEBUG: Insert grammage here
        }

        track.push_back(track.back());
        sector_changed = AdvanceParticle(track.back(), energy_at_next_interaction, distance_to_next_interaction, rnd, current_sector); // DEBUG: We know grammage, so we will pass grammage here

        if (!sector_changed) {
            switch (next_interaction_type) {
                case Stochastic: {
                    track.push_back(track.back());
                    continue_propagation = DoStochasticInteraction(track.back(), get<UTILITY>(current_sector), rnd);
                    break;
                }
                case Decay:
                    track.back().SetType(InteractionType::Decay);
                    continue_propagation = false;
                    break;
                case MaxDistance:
                case MinimalE: {
                    continue_propagation = false;
                    break;
                }
                case ApproachingSector: {
                    continue_propagation = true; //?
                    break;
                }
            }
        } else {
            track.push_back(track.back());
            current_sector = ChooseCurrentSector(track.back().GetPosition(), track.back().GetDirection());
        }

        if(track.back().GetEnergy() <= min_energy)
            continue_propagation = false;
    }

    return track;
}

// Private methods

bool Propagator::DoStochasticInteraction(DynamicData& p_cond, PropagationUtility& utility, std::function<double()> rnd) {
    auto stochastic_cross = utility.TypeInteraction(p_cond.GetEnergy(), rnd);
    double loss_energy = utility.EnergyStochasticloss(*stochastic_cross, p_cond.GetEnergy(), rnd);

    auto deflection_angles = stochastic_cross->StochasticDeflection(p_cond.GetEnergy(), loss_energy); //TODO: ugly
    p_cond.DeflectDirection(get<0>(deflection_angles), get<1>(deflection_angles));

    p_cond.SetEnergy(p_cond.GetEnergy() - loss_energy);
    p_cond.SetType(stochastic_cross->GetParametrization().GetInteractionType());

    return true; //TODO: Add condition for fatal interaction

}

bool Propagator::AdvanceParticle(DynamicData& p_cond, double E_f, double advance_distance, std::function<double()> rnd,
                                 Sector& sector ) {
    assert(advance_distance > 0); // DEBUG: advance_distance is grammage now!
    Vector3D new_position;
    Vector3D new_direction;
    Vector3D mean_direction;
    bool sector_changed = false;

    //TODO: For NoScattering, these random numbers are not used
    const std::array<double, 4> rnd_scattering = {rnd(), rnd(), rnd(), rnd()};

    std::tie(mean_direction, new_direction) = get<UTILITY>(sector).DirectionsScatter(
            advance_distance, p_cond.GetEnergy(), E_f, p_cond.GetDirection(), rnd_scattering); // DEBUG: advance_distance is grammage now (which is ok because Scatter does want grammage now)

    new_position = p_cond.GetPosition() + advance_distance * mean_direction; // DEBUG: Beforehand, we have to convert advance_distance (which is grammage now) to an actual distance using "Correct"

    //TODO: We could also compare the utilities here, but this may lead to errors for weird combined geometries
    //TODO: For hollow geometries, this can be problematic. We probably have to remove them.
    //TODO: Catch behaviour for no defined sector
    if (get<GEOMETRY>(sector) != get<GEOMETRY>(ChooseCurrentSector(new_position, new_direction))) {
        sector_changed = true;
        double old_distance_to_border; // DEBUG this is still a distance
        std::cout << "HIER" << std::endl;
        do {
            old_distance_to_border = advance_distance; //DEBUG: This is a distance now
            advance_distance = CalculateDistanceToBorder(p_cond.GetPosition(), mean_direction, *get<GEOMETRY>(sector)); //DEBUG: Distance
            std::cout << "H: " << old_distance_to_border << ", " << advance_distance << std::endl;
            E_f = get<UTILITY>(sector).EnergyDistance(p_cond.GetEnergy(), advance_distance); // DEBUG: We need to turn the distance to grammage first (using integrate) so we can use EnergyDistance (which expects grammage)
            std::tie(mean_direction, new_direction) = get<UTILITY>(sector).DirectionsScatter(
                    advance_distance, p_cond.GetEnergy(), E_f, p_cond.GetDirection(), rnd_scattering); // DEBUG: This expects the grammage we calculated above
            std::cout << std::setprecision(15) << "Step: " << old_distance_to_border << " vs " << advance_distance << std::endl;
        } while (std::abs(old_distance_to_border - advance_distance) > PARTICLE_POSITION_RESOLUTION); // DEBUG: This is both still distance
        new_position = p_cond.GetPosition() + advance_distance * mean_direction; // DEBUG: This is still a distance
    }

    p_cond.SetPosition(new_position);
    p_cond.SetDirection(new_direction);
    p_cond.SetPropagatedDistance(p_cond.GetPropagatedDistance() + advance_distance); // DEBUG: This is still a distance
    p_cond.SetTime(p_cond.GetTime() + get<UTILITY>(sector).TimeElapsed(p_cond.GetEnergy(), E_f, advance_distance)); // DEBUG: advance_distance is a distance; we need to correct the time for density
    p_cond.SetEnergy(get<UTILITY>(sector).EnergyRandomize(p_cond.GetEnergy(), E_f, rnd));
    p_cond.SetType(InteractionType::ContinuousEnergyLoss);

    return sector_changed;
}

double Propagator::CalculateDistanceToBorder(const Vector3D& position, const Vector3D& direction, const Geometry& current_geometry) {
    double distance_to_border = 0;
    double tmp_distance_to_border = 0;

    distance_to_border = current_geometry.DistanceToBorder(position, direction).first;
    //TODO: These lines should be redundant as soon as the adaptive step size is implemented
    for(auto i = 0; i < sector_list.size(); i++) {
        if(get<GEOMETRY>(sector_list[i])->GetHierarchy() > current_geometry.GetHierarchy()) {
            tmp_distance_to_border = get<GEOMETRY>(sector_list[i])->DistanceToBorder(position, direction).first;
            distance_to_border = std::min(distance_to_border, tmp_distance_to_border);
        }
    }

    return distance_to_border;
}

int Propagator::maximize(const std::array<double, 3>& InteractionEnergies) {
    for (const auto& E : InteractionEnergies) //DEBUG
        std::cout << E << ", ";
    std::cout << std::endl;
    auto max_element_ref = std::max_element(InteractionEnergies.begin(), InteractionEnergies.end());
    return std::distance(InteractionEnergies.begin(), max_element_ref);
}

Sector Propagator::ChooseCurrentSector(const Vector3D& particle_position, const Vector3D& particle_direction) {
    std::vector<int> current_sector_indices;
    for(auto i = 0; i < sector_list.size(); i++) {
        if (get<GEOMETRY>(sector_list[i])->IsInside(particle_position, particle_direction)) {
            current_sector_indices.push_back(i);
        }
    }

    if(current_sector_indices.size() == 0) {
        throw std::logic_error("Propagator: No sector defined at current particle position.");
    } else if (current_sector_indices.size() == 1) {
        return sector_list[current_sector_indices[0]];
    }

    double current_index = current_sector_indices.back();
    for (auto i: current_sector_indices) {
        if (get<GEOMETRY>(sector_list[current_index])->GetHierarchy() <
            get<GEOMETRY>(sector_list[i])->GetHierarchy()) {
            current_index = i;
        }
    }

    return sector_list.at(current_index);
}

// Init methods

nlohmann::json Propagator::ParseConfig(const std::string& config_file){
    nlohmann::json json_config;
    try {
        std::string expanded_config_file_path = Helper::ResolvePath(config_file, true);
        std::ifstream infilestream(expanded_config_file_path);
        infilestream >> json_config;
    } catch (const nlohmann::json::parse_error& e) {
        log_fatal("Unable parse \"%s\" as json file", config_file.c_str());
    }
    return json_config;
};

void Propagator::InitializeSectorFromJSON(const ParticleDef& p_def, const nlohmann::json& json_sector, GlobalSettings global){

    bool do_interpolation = json_sector.value("do_interpolation", global.do_interpolation);
    bool do_exact_time = json_sector.value("exact_time", global.do_exact_time);
    std::string scattering = json_sector.value("scattering", global.scattering);

    std::shared_ptr<Medium> medium = global.medium;
    if (json_sector.contains("medium")) {
        medium = CreateMedium(json_sector["medium"]);
    } else if (medium == nullptr) {
        throw std::invalid_argument("Neither a specific Sector medium nor a global medium is defined.");
    }

    std::shared_ptr<EnergyCutSettings> cuts = global.cuts;
    if (json_sector.contains("cuts")) {
        cuts = std::make_shared<EnergyCutSettings>(EnergyCutSettings(json_sector["cuts"]));
    } else if (cuts == nullptr) {
        throw std::invalid_argument("Neither a specific Sector EnergyCut nor a global EnergyCut is defined.");
    }

    CrossSectionList crosssections;
    if (json_sector.contains("CrossSections")) {
        crosssections = CreateCrossSectionList(json_sector["CrossSections"], medium, cuts, do_interpolation);
    } else if (global.cross != "") {
        crosssections = CreateCrossSectionList(global.cross, medium, cuts, do_interpolation);
    } else {
        crosssections = GetStdCrossSections(medium, cuts, p_def);
    }

    // Create Geometry and add them to sector_list with corresponding geometries
    auto utility = CreateUtility(crosssections, medium, cuts->GetContRand(), do_interpolation, do_exact_time, scattering);
    if (json_sector.contains("geometries")) {
        assert(json_sector["geometries"].is_array());

        nlohmann::json density_distr_default = {{"density_distr_type", "homogeneous"}}; //if no density_distr defined, fallback to homogeneous medium
        auto density_distr_sector = json_sector.value("density_distr", density_distr_default);

        for (const auto& json_geometry : json_sector.at("geometries")) {
            auto density_distr_json = json_geometry.value("density_distr", density_distr_sector);
            sector_list.push_back(std::make_tuple(
                    CreateGeometry(json_geometry), utility, CreateDensityDistribution(density_distr_json)));
        }
    } else {
        throw std::invalid_argument("At least one geometry must be defined for each sector");
    }
}

PropagationUtility Propagator::CreateUtility(CrossSectionList crosssections, std::shared_ptr<Medium> medium, bool do_cont_rand,
                                             bool do_interpolation, bool do_exact_time, std::string scattering ) {
    PropagationUtility::Collection def;
    // Interaction
    if (do_interpolation) {
        def.interaction_calc = std::make_shared<InteractionBuilder<UtilityInterpolant>>(crosssections);
    } else {
        def.interaction_calc = std::make_shared<InteractionBuilder<UtilityIntegral>>(crosssections);
    }

    // Displacement
    if (do_interpolation) {
        def.displacement_calc = std::make_shared<DisplacementBuilder<UtilityInterpolant>>(crosssections);
    } else {
        def.displacement_calc = std::make_shared<DisplacementBuilder<UtilityIntegral>>(crosssections);
    }

    // Time
    if (do_exact_time) {
        if (do_interpolation) {
            def.time_calc = std::make_shared<ExactTimeBuilder<UtilityInterpolant>>(crosssections, p_def);
        } else {
            def.time_calc = std::make_shared<ExactTimeBuilder<UtilityIntegral>>(crosssections, p_def);
        }
    } else {
        def.time_calc = std::make_shared<ApproximateTimeBuilder>();
    }

    // Scattering
    if (scattering != "" and scattering != "NoScattering") {
        if (do_interpolation) {
            //TODO: Making a shared_ptr out of a raw pointer is bad-practise (and looks ugly), the underlying factory should be changed (jm)
            auto aux = ScatteringFactory::Get().CreateScattering(
                    scattering, p_def, medium, interpol_def_global,
                    std::unique_ptr<CrossSectionList>(new CrossSectionList(crosssections)));
            def.scattering = std::shared_ptr<Scattering>(aux);
        } else {
            auto aux = ScatteringFactory::Get().CreateScattering(
                    scattering, p_def, medium, std::unique_ptr<CrossSectionList>(new CrossSectionList(crosssections)));
            def.scattering = std::shared_ptr<Scattering>(aux);
        }
    }

    //Decay
    if (isfinite(p_def.lifetime)) {
        if (do_interpolation) {
            def.decay_calc = std::make_shared<DecayBuilder<UtilityInterpolant>>(crosssections, p_def);
        }
        else {
            def.decay_calc = std::make_shared<DecayBuilder<UtilityIntegral>>(crosssections, p_def);
        }
    }

    //ContRand
    if (do_cont_rand) {
        if (do_interpolation) {
            def.cont_rand = std::make_shared<ContRandBuilder<UtilityInterpolant>>(crosssections);
        } else {
            def.cont_rand = std::make_shared<ContRandBuilder<UtilityIntegral>>(crosssections);
        }
    }

    return def;
}

CrossSectionList Propagator::CreateCrossSectionList(const nlohmann::json& config_cross, std::shared_ptr<Medium> medium,
                                                    std::shared_ptr<EnergyCutSettings> cuts, bool do_interpolation) {
    CrossSectionList crosssections;
    std::shared_ptr<InterpolationDef> interpol_def = nullptr;
    if( do_interpolation ) {
        interpol_def = interpol_def_global;
    }

    //Annihilation
    if (config_cross.contains("annihilation")) {
        AnnihilationFactory::Definition annihilation_def(config_cross["annihilation"]);
        crosssections.emplace_back(
                AnnihilationFactory::Get().CreateAnnihilation(p_def, medium, annihilation_def, interpol_def));
    }

    //Bremsstrahlung
    if (config_cross.contains("brems")) {
        BremsstrahlungFactory::Definition brems_def(config_cross["brems"]);
        crosssections.emplace_back(
                BremsstrahlungFactory::Get().CreateBremsstrahlung(p_def, medium, cuts, brems_def, interpol_def));
    }

    //Compton
    if (config_cross.contains("compton")) {
        ComptonFactory::Definition compton_def(config_cross["compton"]);
        crosssections.emplace_back(
                ComptonFactory::Get().CreateCompton(p_def, medium, cuts, compton_def, interpol_def));
    }

    //EpairProduction
    if (config_cross.contains("epair")) {
        EpairProductionFactory::Definition epair_def(config_cross["epair"]);
        crosssections.emplace_back(
                EpairProductionFactory::Get().CreateEpairProduction(p_def, medium, cuts, epair_def, interpol_def));
    }

    //Ionization
    if (config_cross.contains("ioniz")) {
        IonizationFactory::Definition ioniz_def(config_cross["ioniz"]);
        crosssections.emplace_back(
                IonizationFactory::Get().CreateIonization(p_def, medium, cuts, ioniz_def, interpol_def));
    }

    //MupairProduction
    if (config_cross.contains("mupair")) {
        MupairProductionFactory::Definition mupair_def(config_cross["mupair"]);
        crosssections.emplace_back(
                MupairProductionFactory::Get().CreateMupairProduction(p_def, medium, cuts, mupair_def, interpol_def));
    }

    //Photonulear
    if (config_cross.contains("photonuclear")) {
        PhotonuclearFactory::Definition photonuclear_def(config_cross["photonuclear"]);
        crosssections.emplace_back(
                PhotonuclearFactory::Get().CreatePhotonuclear(p_def, medium, cuts, photonuclear_def, interpol_def));
    }

    //PhotoPair
    if (config_cross.contains("photopair")) {
        PhotoPairFactory::Definition photopair_def(config_cross["photopair"]);
        crosssections.emplace_back(
                PhotoPairFactory::Get().CreatePhotoPair(p_def, medium, photopair_def, interpol_def));
    }

    //WeakInteraction
    if (config_cross.contains("weak")) {
        WeakInteractionFactory::Definition weak_def(config_cross["weak"]);
        crosssections.emplace_back(
                WeakInteractionFactory::Get().CreateWeakInteraction(p_def, medium, weak_def, interpol_def));
    }

    if(crosssections.size()==0){
        throw std::invalid_argument("No crosssections could be initialized for sector");
    }
    return crosssections;
}

Propagator::GlobalSettings::GlobalSettings(const nlohmann::json& config_global) {
    if (config_global.contains("medium")) {
        medium = CreateMedium(config_global);
    }

    if (config_global.contains("CrossSections")) {
        cross = config_global["CrossSections"];
    } else {
        cross = "";
    }

    if (config_global.contains("cuts")) {
        cuts = std::make_shared<EnergyCutSettings>(config_global["cuts"]);
    }

    //Read global propagation settings
    do_exact_time = config_global.value("exact_time", true);
    do_interpolation = config_global.value("do_interpolation", true);
    scattering = config_global.value("scattering", "");
}

Propagator::GlobalSettings::GlobalSettings() {
    medium = nullptr;
    cross = nullptr;
    do_exact_time = true;
    do_interpolation = true;
    scattering = "";
}
