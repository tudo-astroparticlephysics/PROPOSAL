/*
 * Propagator.cxx
 *
 *  Created on: 23.04.2013
 *      Author: koehne
 */

// #include <cmath>

#include <PROPOSAL/crossection/factories/PhotoPairFactory.h>
#include <fstream>
#include <memory>

#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/geometry/Box.h"
#include "PROPOSAL/geometry/Cylinder.h"
#include "PROPOSAL/geometry/GeometryFactory.h"
#include "PROPOSAL/geometry/Sphere.h"

#include "PROPOSAL/medium/MediumFactory.h"

#include "PROPOSAL/particle/Particle.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/json.hpp"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/math/RandomGenerator.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Global defaults
// ------------------------------------------------------------------------- //

const int Propagator::global_seed_ = 0;
const double Propagator::global_ecut_inside_ = 500;
const double Propagator::global_ecut_infront_ = -1;
const double Propagator::global_ecut_behind_ = -1;
const double Propagator::global_vcut_inside_ = -1;
const double Propagator::global_vcut_infront_ = 0.001;
const double Propagator::global_vcut_behind_ = -1;
const double Propagator::global_cont_inside_ = false;
const double Propagator::global_cont_infront_ = true;
const double Propagator::global_cont_behind_ = false;
const bool Propagator::do_interpolation_ = true;
const bool Propagator::uniform_ = true;

// ------------------------------------------------------------------------- //
// Constructors & destructor
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
Propagator::Propagator(
    const std::vector<Sector*>& sectors, const Geometry& geometry) try
    : current_sector_(NULL),
      particle_def_(sectors.at(0)->GetParticleDef()),
      detector_(geometry.clone()),
      /* entry_condition_(DynamicData()), */
      /* exit_condition_(DynamicData()), */
      closest_approach_condition_(DynamicData()),
      elost_(0.) {
    // --------------------------------------------------------------------- //
    // Check if all ParticleDefs are the same
    // --------------------------------------------------------------------- //

    for (auto sector : sectors) {
        if (sector->GetParticleDef() != particle_def_) {
            log_fatal("The particle definitions of the sectors must be "
                      "identical for proper propagation!");
        } else {
            sectors_.push_back(new Sector(*sector));
        }
    }

    current_sector_ = sectors_.at(0);
} catch (const std::out_of_range& ex) {
    log_fatal("No Sectors are provided for the Propagator!");
}

// ------------------------------------------------------------------------- //
Propagator::Propagator(const ParticleDef& particle_def,
    const std::vector<Sector::Definition>& sector_defs,
    const Geometry& geometry)
    : particle_def_(particle_def)
    , detector_(geometry.clone())
    /* , entry_condition_(DynamicData()) */
    /* , exit_condition_(DynamicData()) */
    , closest_approach_condition_(DynamicData())
    , elost_(0.)
{
    for (auto def : sector_defs) {
        sectors_.push_back(new Sector(particle_def, def));
    }

    try {
        current_sector_ = sectors_.at(0);
    } catch (const std::out_of_range& ex) {
        log_fatal("No Sectors are provided for the Propagator!");
    }
}

// ------------------------------------------------------------------------- //
Propagator::Propagator(const ParticleDef& particle_def,
    const std::vector<Sector::Definition>& sector_defs,
    const Geometry& geometry, const InterpolationDef& interpolation_def)
    : particle_def_(particle_def)
    , detector_(geometry.clone())
    /* , entry_condition_(DynamicData()) */
    /* , exit_condition_(DynamicData()) */
    , closest_approach_condition_(DynamicData())
    , elost_(0.)
{
    for (auto def : sector_defs) {
        sectors_.push_back(new Sector(particle_def, def, interpolation_def));
    }

    try {
        current_sector_ = sectors_.at(0);
    } catch (const std::out_of_range& ex) {
        log_fatal("No Sectors are provided for the Propagator!");
    }
}

// ------------------------------------------------------------------------- //
Propagator::Propagator(const Propagator& propagator)
    : sectors_(propagator.sectors_.size(), NULL)
    , current_sector_(NULL)
    , particle_def_(propagator.particle_def_)
    , detector_(propagator.detector_->clone())
    /* , entry_condition_(propagator.entry_condition_) */
    /* , exit_condition_(propagator.exit_condition_) */
    , closest_approach_condition_(propagator.closest_approach_condition_)
    , elost_(propagator.elost_)
{
    for (unsigned int i = 0; i < propagator.sectors_.size(); ++i) {
        sectors_[i] = new Sector(particle_def_, *propagator.sectors_[i]);

        if (propagator.sectors_[i] == propagator.current_sector_) {
            current_sector_ = sectors_[i];
        }
    }
}

// ------------------------------------------------------------------------- //
Propagator::Propagator(
    const ParticleDef& particle_def, const std::string& config_file)
    : current_sector_(NULL)
    , particle_def_(particle_def)
    , detector_(NULL)
    /* , entry_condition_(DynamicData()) */
    /* , exit_condition_(DynamicData()) */
    , closest_approach_condition_(DynamicData())
    , elost_(0.)
{
    int global_seed = global_seed_;
    bool do_interpolation = do_interpolation_;
    bool uniform = uniform_;

    InterpolationDef interpolation_def;
    Sector::Definition sec_def_global;

    // Create the json parser
    nlohmann::json json_config;
    try {
        std::string expanded_config_file_path
            = Helper::ResolvePath(config_file, true);
        std::ifstream infilestream(expanded_config_file_path);
        infilestream >> json_config;
    } catch (const nlohmann::json::parse_error& e) {
        log_fatal("Unable parse \"%s\" as json file", config_file.c_str());
    }

    // set global settings
    if (json_config.find("global") == json_config.end()) {
        log_fatal("The 'globals' option is not set. Use default options.");
    }
    if (!json_config["global"].is_object()) {
        log_fatal("Invalid input for option 'nodes_propagate'. Expected a json "
                  "object.");
    }
    nlohmann::json json_global = json_config["global"];
    std::string json_global_str = json_global.dump();

    // set the seed of the random number generator
    if (json_global.find("seed") != json_global.end()) {
        if (json_global["seed"].is_number()) {
            global_seed = json_global["seed"].get<int>();
            // this can throw the error (nlohmann::detail::type_error& e) if no
            // numbercheck
        } else {
            log_fatal("Invalid input for option 'seed'. Expected a number.");
        }
    }

    RandomGenerator::Get().SetSeed(global_seed);
    log_info("Seed of the default random generator set to %i", global_seed);

    // Read in global cut and continuous randomization options
    std::string cut_object_str;
    nlohmann::json cuts_infront_object;
    cut_object_str = ParseCutSettings(json_global_str, "cuts_infront",
        global_ecut_infront_, global_vcut_infront_, global_cont_infront_);
    cuts_infront_object = nlohmann::json::parse(cut_object_str);

    nlohmann::json cuts_inside_object;
    cut_object_str = ParseCutSettings(json_global_str, "cuts_inside",
        global_ecut_inside_, global_vcut_inside_, global_cont_inside_);
    cuts_inside_object = nlohmann::json::parse(cut_object_str);

    nlohmann::json cuts_behind_object;
    cut_object_str = ParseCutSettings(json_global_str, "cuts_behind",
        global_ecut_behind_, global_vcut_behind_, global_cont_behind_);
    cuts_behind_object = nlohmann::json::parse(cut_object_str);

    // Read in interpolation options
    if (json_global.find("interpolation") != json_global.end()) {
        if (json_global["interpolation"].is_object()) {
            if (json_global["interpolation"].find("do_interpolation")
                != json_global["interpolation"].end()) {
                if (json_global["interpolation"]["do_interpolation"]
                        .is_boolean()) {
                    do_interpolation
                        = json_global["interpolation"]["do_interpolation"];
                } else {
                    log_fatal("Invalid input for option 'do_interpolation'. "
                              "Expected a bool.");
                }
            } else {
                log_debug("The 'do_interpolation' option is not set. Use "
                          "default (true)");
            }
            if (do_interpolation) {
                std::string json_str = json_global["interpolation"].dump();
                interpolation_def = CreateInterpolationDef(json_str);
            } else {
                log_info("Integration instead of interpolation is chosen. The "
                         "propagation is now extremely slow.");
            }
        } else {
            log_fatal("Invalid input for global option 'interpolation'. "
                      "Expected a json object.");
        }
    } else {
        log_debug("The 'interpolation' option is not set. Use default "
                  "interpolation settings.");
    }

    // Read in global sector definition
    sec_def_global = CreateSectorDefinition(json_global_str);

    // read in the uniform flag
    if (json_global.find("uniform") != json_global.end()) {
        if (json_global["uniform"].is_boolean()) {
            uniform = json_global["uniform"].get<bool>();
        } else {
            log_fatal("Invalid input for option 'uniform' (uniform phasespace "
                      "sampling for decays). Expected a bool.");
        }
    } else {
        log_debug("The 'uniform' option (uniform phasespace sampling for "
                  "decays) is not set. Use default (true)");
    }
    particle_def_.decay_table.SetUniformSampling(uniform);

    // Parse detector geometry
    if (json_config.find("detector") != json_config.end()) {
        std::string json_str = json_config["detector"].dump();
        detector_ = ParseGeometryConifg(json_str);
    } else {
        log_fatal("You need to specify a detector geometry.");
    }

    // Read in all sector definitions
    nlohmann::json json_sectors;

    if (json_config.find("sectors") != json_config.end()) {
        if (json_config["sectors"].is_array()) {
            if (1 > json_config["sectors"].size()) {
                log_fatal("You need to specify at least one sector.");
            }
        } else {
            log_fatal("Invalid input for option 'sectors'. Expected an array "
                      "of json objects.");
        }
    } else {
        log_fatal("The 'sectors' option is not set.");
    }

    for (size_t idx = 0; idx < json_config["sectors"].size(); idx++) {
        json_sectors = json_config["sectors"][idx];
        if (!json_sectors.is_object()) {
            log_fatal("Invalid input for an object in 'sectors'. Each sector "
                      "must be a json object.");
        }
        std::string json_sector_str = json_sectors.dump();

        // Create medium
        MediumFactory::Definition medium_def;
        if (json_sectors.find("density_correction") != json_sectors.end()) {
            if (json_sectors["density_correction"].is_number()) {
                medium_def.density_correction
                    = json_sectors["density_correction"].get<double>();
            } else {
                log_fatal("Invalid input for option 'density_correction'. "
                          "Expected a number.");
            }
        } else {
            log_debug("The 'density_correction' option is not set. Use default "
                      "(1.0)");
        }

        std::string medium_name
            = MediumFactory::Get().GetStringFromEnum(medium_def.type);
        if (json_sectors.find("medium") != json_sectors.end()) {
            if (json_sectors["medium"].is_string()) {
                medium_name = json_sectors["medium"].get<std::string>();
            } else {
                log_fatal(
                    "Invalid input for option 'medium'. Expected a string.");
            }
        } else {
            log_debug("The 'medium' option is not set. Use default (Water)");
        }
        Medium* med = MediumFactory::Get().CreateMedium(
            medium_name, medium_def.density_correction);

        // Create Geometry
        Geometry* geometry = NULL;

        if (json_sectors.find("geometry") != json_sectors.end()) {
            std::string json_str = json_sectors["geometry"].dump();
            geometry = ParseGeometryConifg(json_str);
        } else {
            log_fatal("You need to specify a geometry for each sector");
        }

        double hierarchy = geometry->GetHierarchy();
        if (json_sectors.find("hierarchy") != json_sectors.end()) {
            if (json_sectors["hierarchy"].is_number()) {
                hierarchy = json_sectors["hierarchy"].get<int>();
            } else {
                log_fatal(
                    "Invalid input for option 'hierarchy'. Expected a number.");
            }
        } else {
            log_debug("The 'hierarchy' option is not set. Use default (0)");
        }
        geometry->SetHierarchy(hierarchy);

        // Use global options in case they will not be overridden
        Sector::Definition sec_def_infront = sec_def_global;
        sec_def_infront.location = Sector::ParticleLocation::InfrontDetector;
        sec_def_infront.SetMedium(*med);
        sec_def_infront.SetGeometry(*geometry);

        Sector::Definition sec_def_inside = sec_def_global;
        sec_def_inside.location = Sector::ParticleLocation::InsideDetector;
        sec_def_inside.SetMedium(*med);
        sec_def_inside.SetGeometry(*geometry);

        Sector::Definition sec_def_behind = sec_def_global;
        sec_def_behind.location = Sector::ParticleLocation::BehindDetector;
        sec_def_behind.SetMedium(*med);
        sec_def_behind.SetGeometry(*geometry);

        nlohmann::json json_cutsettings;

        // cut settings infront
        cut_object_str = ParseCutSettings(json_sector_str, "cuts_infront",
            cuts_infront_object["e_cut"].get<double>(),
            cuts_infront_object["v_cut"].get<double>(),
            cuts_infront_object["cont_rand"].get<bool>());

        json_cutsettings = nlohmann::json::parse(cut_object_str);
        sec_def_infront.cut_settings.SetEcut(
            json_cutsettings["e_cut"].get<double>());
        sec_def_infront.cut_settings.SetVcut(
            json_cutsettings["v_cut"].get<double>());
        sec_def_infront.do_continuous_randomization
            = json_cutsettings["cont_rand"].get<bool>();

        // cut settings inside
        cut_object_str = ParseCutSettings(json_sector_str, "cuts_inside",
            cuts_inside_object["e_cut"].get<double>(),
            cuts_inside_object["v_cut"].get<double>(),
            cuts_inside_object["cont_rand"].get<bool>());

        json_cutsettings = nlohmann::json::parse(cut_object_str);
        sec_def_inside.cut_settings.SetEcut(
            json_cutsettings["e_cut"].get<double>());
        sec_def_inside.cut_settings.SetVcut(
            json_cutsettings["v_cut"].get<double>());
        sec_def_inside.do_continuous_randomization
            = json_cutsettings["cont_rand"].get<bool>();

        // cut settings behind
        cut_object_str = ParseCutSettings(json_sector_str, "cuts_behind",
            cuts_behind_object["e_cut"].get<double>(),
            cuts_behind_object["v_cut"].get<double>(),
            cuts_behind_object["cont_rand"].get<bool>());

        json_cutsettings = nlohmann::json::parse(cut_object_str);
        sec_def_behind.cut_settings.SetEcut(
            json_cutsettings["e_cut"].get<double>());
        sec_def_behind.cut_settings.SetVcut(
            json_cutsettings["v_cut"].get<double>());
        sec_def_behind.do_continuous_randomization
            = json_cutsettings["cont_rand"].get<bool>();

        if (do_interpolation) {
            sectors_.push_back(
                new Sector(particle_def_, sec_def_infront, interpolation_def));
            sectors_.push_back(
                new Sector(particle_def_, sec_def_inside, interpolation_def));
            sectors_.push_back(
                new Sector(particle_def_, sec_def_behind, interpolation_def));
        } else {
            sectors_.push_back(new Sector(particle_def_, sec_def_infront));
            sectors_.push_back(new Sector(particle_def_, sec_def_inside));
            sectors_.push_back(new Sector(particle_def_, sec_def_behind));
        }

        delete geometry;
        delete med;
    }
}

Propagator::~Propagator()
{
    for (auto sector : sectors_) {
        delete sector;
    }

    sectors_.clear();

    delete detector_;
}

// ------------------------------------------------------------------------- //
// Operators
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
bool Propagator::operator==(const Propagator& propagator) const
{
    if (*detector_ != *propagator.detector_) {
        return false;
    }
    if (sectors_.size() != propagator.sectors_.size()) {
        return false;
    }
    /* if (entry_condition_ != propagator.entry_condition_) { */
    /*     return false; */
    /* } */
    /* if (exit_condition_ != propagator.exit_condition_) { */
    /*     return false; */
    /* } */
    if (closest_approach_condition_ != propagator.closest_approach_condition_) {
        return false;
    }
    if (elost_ != propagator.elost_) {
        return false;
    } else {
        for (unsigned int i = 0; i < sectors_.size(); ++i) {
            if (*sectors_[i] != *propagator.sectors_[i]) {
                return false;
            }
        }
    }

    return true;
}

// ------------------------------------------------------------------------- //
bool Propagator::operator!=(const Propagator& propagator) const
{
    return !(*this == propagator);
}

// ------------------------------------------------------------------------- //
// Public member functions
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
Secondaries Propagator::Propagate(
    DynamicData initial_condition, double max_distance, double minimal_energy)
{
    double distance = 0;
    double distance_to_closest_approach = 0;

    Secondaries secondaries_(std::make_shared<ParticleDef>(particle_def_));
    /* secondaries_.reserve(static_cast<size_t>(produced_particle_moments_.first
     */
    /*     + 2 * std::sqrt(produced_particle_moments_.second))); */

    // These two variables are needed to calculate the energy loss inside the
    // detector energy_at_entry_point is initialized with the current energy
    // because this is a reasonable value for particle which starts inside the
    // detector They should be set below, so there is no need to init them here.
    // But for safetiness, if an edge case is not considered
    // one could include them.
    // particle_.SetEntryEnergy(particle_.GetEnergy());
    // particle_.SetExitEnergy(particle_.GetMass());

    // TODO: what to do with the entry, exit closest approach point?
    // DynamicData entry_condition;
    // DynamicData exit_condition;
    // DynamicData closest_approach_condition;

    bool starts_in_detector = detector_->IsInside(
        initial_condition.GetPosition(), initial_condition.GetDirection());
    if (starts_in_detector) {
        secondaries_.SetEntryPoint(initial_condition);
        distance_to_closest_approach = detector_->DistanceToClosestApproach(
            initial_condition.GetPosition(), initial_condition.GetDirection());
        if (distance_to_closest_approach < 0) {
            closest_approach_condition_.SetPosition(
                initial_condition.GetPosition());
            closest_approach_condition_.SetEnergy(
                initial_condition.GetEnergy());
            closest_approach_condition_.SetTime(initial_condition.GetTime());
        }
    }
    bool is_in_detector = false;
    bool was_in_detector = false;
    bool propagationstep_till_closest_approach = false;
    bool already_reached_closest_approach = false;

    std::unique_ptr<DynamicData> p_condition(
        new DynamicData(initial_condition));
    while (1) {
        ChooseCurrentSector(
            p_condition->GetPosition(), p_condition->GetDirection());

        if (current_sector_ == nullptr) {
            log_info("particle reached the border");
            break;
        }

        // Check if have to propagate the particle_ through the whole sector
        // or only to the sector border
        distance = CalculateEffectiveDistance(
            p_condition->GetPosition(), p_condition->GetDirection());

        if (already_reached_closest_approach == false) {
            distance_to_closest_approach = detector_->DistanceToClosestApproach(
                p_condition->GetPosition(), p_condition->GetDirection());
            if (distance_to_closest_approach > 0) {
                if (distance_to_closest_approach < distance) {
                    already_reached_closest_approach = true;

                    if (std::abs(distance_to_closest_approach)
                        < GEOMETRY_PRECISION) {
                        closest_approach_condition_.SetPosition(
                            p_condition->GetPosition());
                        closest_approach_condition_.SetEnergy(
                            p_condition->GetEnergy());
                        closest_approach_condition_.SetTime(
                            p_condition->GetTime());
                    } else {
                        distance = distance_to_closest_approach;
                        propagationstep_till_closest_approach = true;
                    }
                }
            }
        }

        is_in_detector = detector_->IsInside(
            p_condition->GetPosition(), p_condition->GetDirection());
        // entry point of the detector
        if (!starts_in_detector && !was_in_detector && is_in_detector) {
            secondaries_.SetEntryPoint(*p_condition);

            was_in_detector = true;
        }
        // exit point of the detector
        else if (was_in_detector && !is_in_detector) {
            secondaries_.SetExitPoint(*p_condition);

            // we don't want to run in this case a second time so we set
            // was_in_detector to false
            was_in_detector = false;

        }
        // if particle_ starts inside the detector we only ant to fill the exit
        // point
        else if (starts_in_detector && !is_in_detector) {
            secondaries_.SetExitPoint(*p_condition);

            // we don't want to run in this case a second time so we set
            // starts_in_detector to false
            starts_in_detector = false;
        }
        if (max_distance <= p_condition->GetPropagatedDistance() + distance) {
            distance = max_distance - p_condition->GetPropagatedDistance();
        }

        Secondaries sector_secondaries = current_sector_->Propagate(
            *p_condition, distance, minimal_energy);
        secondaries_.append(sector_secondaries);

        // TODO: this is not god, because the seconday can have other conditions
        p_condition.reset(
            new DynamicData(sector_secondaries.GetSecondaries().back()));

        if (propagationstep_till_closest_approach) {
            closest_approach_condition_.SetPosition(p_condition->GetPosition());
            closest_approach_condition_.SetEnergy(p_condition->GetEnergy());
            closest_approach_condition_.SetTime(p_condition->GetTime());

            propagationstep_till_closest_approach = false;
        }

        if (max_distance <= p_condition->GetPropagatedDistance()
            || p_condition->GetEnergy() <= minimal_energy
            || p_condition->GetTypeId()
                == static_cast<int>(InteractionType::Decay))
            break;
    }
    if (detector_->IsInside(
            p_condition->GetPosition(), p_condition->GetDirection())) {
        secondaries_.SetExitPoint(*p_condition);
    }

    // TODO: what is now with ELost?
    double elost_ = secondaries_.GetExitPoint().GetEnergy()
        - secondaries_.GetExitPoint().GetEnergy();

    secondaries_.DoDecay();

    n_th_call_ += 1.;
    double produced_particles_
        = static_cast<double>(secondaries_.GetNumberOfParticles());
    produced_particle_moments_ = welfords_online_algorithm(produced_particles_,
        n_th_call_, produced_particle_moments_.first,
        produced_particle_moments_.second);

    return secondaries_;
}

// ------------------------------------------------------------------------- //
void Propagator::ChooseCurrentSector(
    const Vector3D& particle_position, const Vector3D& particle_direction)
{
    std::vector<int> crossed_sector;
    crossed_sector.resize(0);

    // Get Location of the detector (Inside/Infront/Behind)
    Geometry::ParticleLocation::Enum detector_location
        = detector_->GetLocation(particle_position, particle_direction);

    for (unsigned int i = 0; i < sectors_.size(); ++i) {
        if (sectors_[i]->GetGeometry()->IsInside(
                particle_position, particle_direction)) {
            if (static_cast<int>(sectors_[i]->GetLocation())
                == static_cast<int>(detector_location))
                crossed_sector.push_back(i);
        }
    }

    // No sector was found
    if (crossed_sector.size() == 0) {
        current_sector_ = nullptr;
        log_warn("There is no sector defined at position [%f, %f, %f] !!!",
            particle_position.GetX(), particle_position.GetY(),
            particle_position.GetZ());
    } else {
        current_sector_ = sectors_[crossed_sector.back()];
    }

    // Choose current sector when multiple sectors are crossed!
    //
    // Choose by hierarchy of Geometry
    // If same hierarchys are available the denser one is chosen
    // If hierarchy and density are the same then the first found is taken.
    //

    for (std::vector<int>::iterator iter = crossed_sector.begin();
         iter != crossed_sector.end(); ++iter) {

        // Current Hierarchy is equal -> Look at the density!
        //
        if (current_sector_->GetGeometry()->GetHierarchy()
            == sectors_[*iter]->GetGeometry()->GetHierarchy()) {
            // Current Density is smaller -> Set the new sector!
            //
            if (current_sector_->GetMedium()->GetCorrectedMassDensity(
                    particle_position)
                < sectors_[*iter]->GetMedium()->GetCorrectedMassDensity(
                      particle_position))
                current_sector_ = sectors_[*iter];
        }

        // Current Hierarchy is smaller -> Set the new sector!
        //
        if (current_sector_->GetGeometry()->GetHierarchy()
            < sectors_[*iter]->GetGeometry()->GetHierarchy())
            current_sector_ = sectors_[*iter];
    }
}

// ------------------------------------------------------------------------- //
double Propagator::CalculateEffectiveDistance(
    const Vector3D& particle_position, const Vector3D& particle_direction)
{
    double distance_to_sector_border = 0;
    double distance_to_detector = 0;

    distance_to_sector_border
        = current_sector_->GetGeometry()
              ->DistanceToBorder(particle_position, particle_direction)
              .first;
    double tmp_distance_to_border;

    Geometry::ParticleLocation::Enum detector_location
        = detector_->GetLocation(particle_position, particle_direction);

    for (std::vector<Sector*>::iterator iter = sectors_.begin();
         iter != sectors_.end(); ++iter) {

        if (static_cast<int>((*iter)->GetLocation())
            == static_cast<int>(detector_location)) {
            if ((*iter)->GetGeometry()->GetHierarchy()
                >= current_sector_->GetGeometry()->GetHierarchy()) {
                tmp_distance_to_border
                    = (*iter)
                          ->GetGeometry()
                          ->DistanceToBorder(
                              particle_position, particle_direction)
                          .first;
                if (tmp_distance_to_border <= 0)
                    continue;
                distance_to_sector_border = std::min(
                    tmp_distance_to_border, distance_to_sector_border);
            }
        }
    }

    distance_to_detector
        = detector_->DistanceToBorder(particle_position, particle_direction)
              .first;

    if (distance_to_detector > 0) {
        return std::min(distance_to_detector, distance_to_sector_border);
    } else {
        return distance_to_sector_border;
    }
}

// ------------------------------------------------------------------------- //
// Private member functions
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
Geometry* Propagator::ParseGeometryConifg(const std::string& json_object_str)
{
    nlohmann::json json_object = nlohmann::json::parse(json_object_str);
    if (!json_object.is_object()) {
        log_fatal("This geometry option should be a json object.");
    }

    std::string origin_str = "origin";
    std::string outer_radius_str = "outer_radius";
    std::string inner_radius_str = "inner_radius";
    std::string length_str = "length";
    std::string width_str = "width";
    std::string height_str = "height";

    std::string warning_str
        = "Geometry %s needs to specify \"%s\" in the config file!";

    double cm_to_meter = 100.0;

    // --------------------------------------------------------------------- //
    // Get geometry from default constructor
    // --------------------------------------------------------------------- //

    Geometry* geometry = NULL;

    if (json_object.find("shape") != json_object.end()) {
        if (json_object["shape"].is_string()) {
            std::string name = json_object["shape"].get<std::string>();
            geometry = GeometryFactory::Get().CreateGeometry(name);
        } else {
            log_fatal("Invalid input for option 'shape'. Expected a string.");
        }
    } else {
        log_fatal(warning_str.c_str(), "", "shape");
    }

    // --------------------------------------------------------------------- //
    // Get the position vector from the property tree
    // --------------------------------------------------------------------- //

    double x = 0;
    double y = 0;
    double z = 0;

    if (json_object.find(origin_str) != json_object.end()) {
        if (json_object[origin_str].is_array()) {
            if (json_object[origin_str].size() == 3) {
                for (size_t idx = 0; idx < json_object[origin_str]
                                               .get<std::vector<double>>()
                                               .size();
                     idx++) {
                    if (json_object[origin_str][idx].is_number()) {
                        double coord
                            = json_object[origin_str][idx].get<double>()
                            * cm_to_meter;

                        switch (idx) {
                        case 0:
                            x = coord;
                            break;
                        case 1:
                            y = coord;
                            break;
                        case 2:
                            z = coord;
                            break;
                        default:
                            // Do nothing
                            break;
                        }
                    } else {
                        log_fatal("Invalid input for option 'origin'. Expected "
                                  "an array consisting of three numbers, but "
                                  "the objects in the array are no numbers.");
                    }
                }
            } else {
                log_fatal("Invalid input for option 'origin'. Expected an "
                          "array consisting of three numbers, but the array "
                          "has not the length three.");
            }
        } else {
            log_fatal(
                "Invalid input for option 'origin'. Expected an array "
                "consisting of three numbers, but the input is not an array.");
        }
    } else {
        log_fatal(warning_str.c_str(), "", "origin");
    }

    Vector3D vec(x, y, z);

    // --------------------------------------------------------------------- //
    // Check type of geometry and set members
    // --------------------------------------------------------------------- //

    if (PROPOSAL::Sphere* sphere = dynamic_cast<PROPOSAL::Sphere*>(geometry)) {
        double radius = 0;
        double inner_radius = 0;

        if (json_object.find(outer_radius_str) != json_object.end()) {
            if (json_object[outer_radius_str].is_number()) {
                radius
                    = json_object[outer_radius_str].get<double>() * cm_to_meter;
            } else {
                log_fatal("Invalid input for option 'outer_radius'. Expected a "
                          "number.");
            }
        } else {
            log_fatal(warning_str.c_str(), sphere->GetName().c_str(),
                outer_radius_str.c_str());
        }

        if (json_object.find(inner_radius_str) != json_object.end()) {
            if (json_object[inner_radius_str].is_number()) {
                inner_radius
                    = json_object[inner_radius_str].get<double>() * cm_to_meter;
            } else {
                log_fatal("Invalid input for option 'inner_radius'. Expected a "
                          "number.");
            }
        } else {
            log_fatal(warning_str.c_str(), sphere->GetName().c_str(),
                inner_radius_str.c_str());
        }

        sphere->SetPosition(vec);
        sphere->SetRadius(radius);
        sphere->SetInnerRadius(inner_radius);

        return sphere;
    } else if (PROPOSAL::Box* box = dynamic_cast<PROPOSAL::Box*>(geometry)) {
        double x = 0;
        double y = 0;
        double z = 0;

        if (json_object.find(length_str) != json_object.end()) {
            if (json_object[length_str].is_number()) {
                x = json_object[length_str].get<double>() * cm_to_meter;
            } else {
                log_fatal(
                    "Invalid input for option 'length'. Expected a number.");
            }
        } else {
            log_fatal(warning_str.c_str(), box->GetName().c_str(),
                length_str.c_str());
        }

        if (json_object.find(width_str) != json_object.end()) {
            if (json_object[width_str].is_number()) {
                y = json_object[width_str].get<double>() * cm_to_meter;
            } else {
                log_fatal(
                    "Invalid input for option 'width'. Expected a number.");
            }
        } else {
            log_fatal(
                warning_str.c_str(), box->GetName().c_str(), width_str.c_str());
        }

        if (json_object.find(height_str) != json_object.end()) {
            if (json_object[height_str].is_number()) {
                z = json_object[height_str].get<double>() * cm_to_meter;
            } else {
                log_fatal(
                    "Invalid input for option 'height'. Expected a number.");
            }
        } else {
            log_fatal(warning_str.c_str(), box->GetName().c_str(),
                height_str.c_str());
        }

        box->SetPosition(vec);
        box->SetX(x);
        box->SetY(y);
        box->SetZ(z);

        return box;
    } else if (PROPOSAL::Cylinder* cylinder
        = dynamic_cast<PROPOSAL::Cylinder*>(geometry)) {
        double radius = 0;
        double inner_radius = 0;
        double z = 0;

        if (json_object.find(outer_radius_str) != json_object.end()) {
            if (json_object[outer_radius_str].is_number()) {
                radius
                    = json_object[outer_radius_str].get<double>() * cm_to_meter;
            } else {
                log_fatal("Invalid input for option 'outer_radius'. Expected a "
                          "number.");
            }
        } else {
            log_fatal(warning_str.c_str(), cylinder->GetName().c_str(),
                outer_radius_str.c_str());
        }

        if (json_object.find(inner_radius_str) != json_object.end()) {
            if (json_object[inner_radius_str].is_number()) {
                inner_radius
                    = json_object[inner_radius_str].get<double>() * cm_to_meter;
            } else {
                log_fatal("Invalid input for option 'inner_radius'. Expected a "
                          "number.");
            }
        } else {
            log_fatal(warning_str.c_str(), cylinder->GetName().c_str(),
                inner_radius_str.c_str());
        }

        if (json_object.find(height_str) != json_object.end()) {
            if (json_object[height_str].is_number()) {
                z = json_object[height_str].get<double>() * cm_to_meter;
            } else {
                log_fatal(
                    "Invalid input for option 'height'. Expected a number.");
            }
        } else {
            log_fatal(warning_str.c_str(), cylinder->GetName().c_str(),
                height_str.c_str());
        }

        cylinder->SetPosition(vec);
        cylinder->SetRadius(radius);
        cylinder->SetInnerRadius(inner_radius);
        cylinder->SetZ(z);

        return cylinder;
    } else {
        log_fatal("Dynamic casts of Geometries failed. Should not end here!");
        return NULL;
    }
}

InterpolationDef Propagator::CreateInterpolationDef(
    const std::string& json_object_str)
{
    // Read in interpolation options
    InterpolationDef interpolation_def;
    nlohmann::json json_object = nlohmann::json::parse(json_object_str);

    if (json_object.find("nodes_propagate") != json_object.end()) {
        if (json_object["nodes_propagate"].is_number()) {
            interpolation_def.nodes_propagate = json_object["nodes_propagate"];
        } else {
            log_fatal("Invalid input for option 'nodes_propagate'. Expected a "
                      "number.");
        }
    } else {
        log_debug(
            "The 'nodes_propagate option' is not set. Use default (1000)");
    }

    if (json_object.find("nodes_continous_randomization")
        != json_object.end()) {
        if (json_object["nodes_continous_randomization"].is_number()) {
            interpolation_def.nodes_continous_randomization
                = json_object["nodes_continous_randomization"];
        } else {
            log_fatal("Invalid input for option "
                      "'nodes_continous_randomization'. Expected a number.");
        }
    } else {
        log_debug("The 'nodes_continous_randomization' option is not set. Use "
                  "default (300)");
    }

    if (json_object.find("nodes_cross_section") != json_object.end()) {
        if (json_object["nodes_cross_section"].is_number()) {
            interpolation_def.nodes_cross_section
                = json_object["nodes_cross_section"];
        } else {
            log_fatal("Invalid input for option 'nodes_cross_section'. "
                      "Expected a number.");
        }
    } else {
        log_debug(
            "The 'nodes_cross_section' option is not set. Use default (100)");
    }

    if (json_object.find("max_node_energy") != json_object.end()) {
        if (json_object["max_node_energy"].is_number()) {
            interpolation_def.max_node_energy = json_object["max_node_energy"];
        } else {
            log_fatal("Invalid input for option 'max_node_energy'. Expected a "
                      "number.");
        }
    } else {
        log_debug(
            "The 'max_node_energy' option is not set. Use default (1e14 MeV)");
    }

    if (json_object.find("do_binary_tables") != json_object.end()) {
        if (json_object["do_binary_tables"].is_boolean()) {
            interpolation_def.do_binary_tables
                = json_object["do_binary_tables"];
        } else {
            log_fatal("Invalid input for option 'do_binary_tables'. Expected a "
                      "bool.");
        }
    } else {
        log_debug(
            "The 'do_binary_tables' option is not set. Use default (true)");
    }

    if (json_object.find("just_use_readonly_path") != json_object.end()) {
        if (json_object["just_use_readonly_path"].is_boolean()) {
            interpolation_def.just_use_readonly_path
                = json_object["just_use_readonly_path"];
        } else {
            log_fatal("Invalid input for option 'just_use_readonly_path'. "
                      "Expected a bool.");
        }
    } else {
        log_debug("The 'just_use_readonly_path' option is not set. Use default "
                  "(false)");
    }

    // Parse to find path to interpolation tables
    if (json_object.find("path_to_tables") != json_object.end()) {
        std::string table_path_str = "";
        if (json_object["path_to_tables"].is_string()) {
            table_path_str = json_object["path_to_tables"];
            interpolation_def.path_to_tables
                = Helper::ResolvePath(table_path_str);
        } else if (json_object["path_to_tables"].is_array()) {
            for (size_t idx = 0; idx < json_object["path_to_tables"]
                                           .get<std::vector<std::string>>()
                                           .size();
                 idx++) {
                if (json_object["path_to_tables"][idx].is_string()) {
                    table_path_str = Helper::ResolvePath(
                        json_object["path_to_tables"][idx].get<std::string>());
                    if (table_path_str != "")
                        break;
                } else {
                    log_fatal("Invalid input for option 'path_to_tables'. "
                              "Expected a string or a list of strings.");
                }
            }
        } else {
            log_fatal("Invalid input for option 'path_to_tables'. Expected a "
                      "string or a list of strings.");
        }

        if (table_path_str != "") {
            interpolation_def.path_to_tables = table_path_str;
            log_info("The writable Path to interpolation tables set to: \"%s\"",
                table_path_str.c_str());
        } else {
            log_warn(
                "No valid writable path to interpolation tables found. Save "
                "tables in memory, if readonly path is also not working!");
        }
    } else {
        log_debug("No writable path to tables set. Use default and save in "
                  "memory, if readonly path is also not working!");
    }

    // Parse to find path to interpolation tables for readonly
    if (json_object.find("path_to_tables_readonly") != json_object.end()) {
        std::string table_path_str = "";
        if (json_object["path_to_tables_readonly"].is_string()) {
            table_path_str = json_object["path_to_tables_readonly"];
            interpolation_def.path_to_tables_readonly
                = Helper::ResolvePath(table_path_str, true);
        } else if (json_object["path_to_tables_readonly"].is_array()) {
            for (size_t idx = 0; idx < json_object["path_to_tables_readonly"]
                                           .get<std::vector<std::string>>()
                                           .size();
                 idx++) {
                if (json_object["path_to_tables_readonly"][idx].is_string()) {
                    table_path_str = Helper::ResolvePath(
                        json_object["path_to_tables_readonly"][idx]
                            .get<std::string>(),
                        true);
                    if (table_path_str != "")
                        break;
                } else {
                    log_fatal(
                        "Invalid input for option 'path_to_tables_readonly'. "
                        "Expected a string or a list of strings.");
                }
            }
        } else {
            log_fatal("Invalid input for option 'path_to_tables_readonly'. "
                      "Expected a string or a list of strings.");
        }

        if (table_path_str != "") {
            interpolation_def.path_to_tables_readonly = table_path_str;
            log_info("Path to readonly interpolation tables set to: \"%s\"",
                table_path_str.c_str());
        } else {
            log_warn("No valid path to readonly interpolation tables found. "
                     "Just looking at writable path_to_tables.");
        }
    } else {
        log_debug("No 'path_to_tables_readonly' option set. Use default and "
                  "look at writable path to tables.");
    }

    return interpolation_def;
}

Sector::Definition Propagator::CreateSectorDefinition(
    const std::string& json_object_str)
{
    Sector::Definition sec_def_global;
    nlohmann::json json_global = nlohmann::json::parse(json_object_str);

    if (json_global.find("brems_multiplier") != json_global.end()) {
        if (json_global["brems_multiplier"].is_number()) {
            sec_def_global.utility_def.brems_def.multiplier
                = json_global["brems_multiplier"].get<double>();
        } else {
            log_fatal("Invalid input for option 'brems_multiplier'. Expected a "
                      "number.");
        }
    } else {
        log_debug("The 'brems_multiplier' option is not set. Use default (%f)",
            sec_def_global.utility_def.brems_def.multiplier);
    }

    if (json_global.find("photo_multiplier") != json_global.end()) {
        if (json_global["photo_multiplier"].is_number()) {
            sec_def_global.utility_def.photo_def.multiplier
                = json_global["photo_multiplier"].get<double>();
        } else {
            log_fatal("Invalid input for option 'photo_multiplier'. Expected a "
                      "number.");
        }

    } else {
        log_debug("The 'photo_multiplier' option is not set. Use default (%f)",
            sec_def_global.utility_def.photo_def.multiplier);
    }

    if (json_global.find("ioniz_multiplier") != json_global.end()) {
        if (json_global["ioniz_multiplier"].is_number()) {
            sec_def_global.utility_def.ioniz_def.multiplier
                = json_global["ioniz_multiplier"].get<double>();
        } else {
            log_fatal("Invalid input for option 'ioniz_multiplier'. Expected a "
                      "number.");
        }
    } else {
        log_debug("The 'ioniz_multiplier' option is not set. Use default (%f)",
            sec_def_global.utility_def.ioniz_def.multiplier);
    }

    if (json_global.find("epair_multiplier") != json_global.end()) {
        if (json_global["epair_multiplier"].is_number()) {
            sec_def_global.utility_def.epair_def.multiplier
                = json_global["epair_multiplier"].get<double>();
        } else {
            log_fatal("Invalid input for option 'epair_multiplier'. Expected a "
                      "number.");
        }
    } else {
        log_debug("The 'epair_multiplier' option is not set. Use default (%f)",
            sec_def_global.utility_def.epair_def.multiplier);
    }

    if (json_global.find("annihilation_multiplier") != json_global.end()) {
        if (json_global["annihilation_multiplier"].is_number()) {
            sec_def_global.utility_def.annihilation_def.multiplier
                = json_global["annihilation_multiplier"].get<double>();
        } else {
            log_fatal("Invalid input for option 'annihilation_multiplier'. "
                      "Expected a number.");
        }
    } else {
        log_debug(
            "The 'annihilation_multiplier' option is not set. Use default (%f)",
            sec_def_global.utility_def.annihilation_def.multiplier);
    }

    if (json_global.find("mupair_multiplier") != json_global.end()) {
        if (json_global["mupair_multiplier"].is_number()) {
            sec_def_global.utility_def.mupair_def.multiplier
                = json_global["mupair_multiplier"].get<double>();
        } else {
            log_fatal("Invalid input for option 'mupair_multiplier'. Expected "
                      "a number.");
        }
    } else {
        log_debug("The 'mupair_multiplier' option is not set. Use default (%f)",
            sec_def_global.utility_def.mupair_def.multiplier);
    }

    if (json_global.find("weak_multiplier") != json_global.end()) {
        if (json_global["weak_multiplier"].is_number()) {
            sec_def_global.utility_def.weak_def.multiplier
                = json_global["weak_multiplier"].get<double>();
        } else {
            log_fatal("Invalid input for option 'weak_multiplier'. Expected a "
                      "number.");
        }
    } else {
        log_debug("The 'weak_multiplier' option is not set. Use default (%f)",
            sec_def_global.utility_def.weak_def.multiplier);
    }

    if (json_global.find("compton_multiplier") != json_global.end()) {
        if (json_global["compton_multiplier"].is_number()) {
            sec_def_global.utility_def.compton_def.multiplier
                = json_global["compton_multiplier"].get<double>();
        } else {
            log_fatal("Invalid input for option 'compton_multiplier'. Expected "
                      "a number.");
        }
    } else {
        log_debug(
            "The 'compton_multiplier' option is not set. Use default (%f)",
            sec_def_global.utility_def.compton_def.multiplier);
    }

    if (json_global.find("photopair_multiplier") != json_global.end()) {
        if (json_global["photopair_multiplier"].is_number()) {
            sec_def_global.utility_def.photopair_def.multiplier
                = json_global["photopair_multiplier"].get<double>();
        } else {
            log_fatal("Invalid input for option 'photopair_multiplier'. "
                      "Expected a number.");
        }
    } else {
        log_debug(
            "The 'photopair_multiplier' option is not set. Use default (%f)",
            sec_def_global.utility_def.photopair_def.multiplier);
    }

    if (json_global.find("lpm") != json_global.end()) {
        if (json_global["lpm"].is_boolean()) {
            sec_def_global.utility_def.epair_def.lpm_effect
                = json_global["lpm"].get<bool>();
            sec_def_global.utility_def.brems_def.lpm_effect
                = json_global["lpm"].get<bool>();
        } else {
            log_fatal("Invalid input for option 'lpm'. Expected a bool.");
        }
    } else {
        log_debug("The 'lpm' option is not set. Use default (true)");
    }

    if (json_global.find("mupair_particle_output") != json_global.end()) {
        if (json_global["mupair_particle_output"].is_boolean()) {
            sec_def_global.utility_def.mupair_def.particle_output
                = json_global["mupair_particle_output"].get<bool>();
        } else {
            log_fatal("Invalid input for option 'mupair_particle_output'. "
                      "Expected a bool.");
        }
    } else {
        log_debug("The 'mupair_particle_output' option is not set. Use default "
                  "(true)");
    }

    if (json_global.find("exact_time") != json_global.end()) {
        if (json_global["exact_time"].is_boolean()) {
            sec_def_global.do_exact_time_calculation
                = json_global["exact_time"].get<bool>();
        } else {
            log_fatal(
                "Invalid input for option 'exact_time'. Expected a bool.");
        }
    } else {
        log_debug("The 'exact_time' option is not set. Use default (true)");
    }

    if (json_global.find("continous_loss_output") != json_global.end()) {
        if (json_global["continous_loss_output"].is_boolean()) {
            sec_def_global.do_continuous_energy_loss_output
                = json_global["continous_loss_output"].get<bool>();
        } else {
            log_fatal("Invalid input for option "
                      "'do_continuous_energy_loss_output'. Expected a bool.");
        }
    } else {
        log_debug("The 'do_continuous_energy_loss_output' option is not set. "
                  "Use default (true)");
    }

    if (json_global.find("stopping_decay") != json_global.end()) {
        if (json_global["stopping_decay"].is_boolean()) {
            sec_def_global.stopping_decay
                = json_global["stopping_decay"].get<bool>();
        } else {
            log_fatal(
                "Invalid input for option 'stopping_decay'. Expected a bool.");
        }
    } else {
        log_debug("The 'stopping_decay' option is not set. Use default (true)");
    }

    if (json_global.find("only_loss_inside_detector") != json_global.end()) {
        if (json_global["only_loss_inside_detector"].is_boolean()) {
            sec_def_global.only_loss_inside_detector
                = json_global["only_loss_inside_detector"].get<bool>();
        } else {
            log_fatal("Invalid input for option 'only_loss_inside_detector'. "
                      "Expected a bool.");
        }
    } else {
        log_debug("The 'only_loss_inside_detector' option is not set. Use "
                  "default (false)");
    }

    if (json_global.find("stochastic_loss_weighting") != json_global.end()) {
        if (json_global["stochastic_loss_weighting"].is_number()) {
            sec_def_global.stochastic_loss_weighting
                = json_global["stochastic_loss_weighting"].get<double>();
            sec_def_global.do_stochastic_loss_weighting
                = 1e-10 < std::abs(sec_def_global.stochastic_loss_weighting);
        } else {
            log_fatal("Invalid input for option 'stochastic_loss_weighting'. "
                      "Expected a number.");
        }
    } else {
        log_debug("The 'stochastic_loss_weighting option' is not set. Use "
                  "default (%f)",
            sec_def_global.stochastic_loss_weighting);
    }

    if (json_global.find("scattering") != json_global.end()) {
        if (json_global["scattering"].is_string()) {
            std::string scattering
                = json_global["scattering"].get<std::string>();
            sec_def_global.scattering_model
                = ScatteringFactory::Get().GetEnumFromString(scattering);
        } else {
            log_fatal(
                "Invalid input for option 'scattering'. Expected a string.");
        }
    } else {
        log_debug("The 'scattering' option is not set. Use default (%s)",
            ScatteringFactory::Get()
                .GetStringFromEnum(sec_def_global.scattering_model)
                .c_str());
    }

    if (json_global.find("brems") != json_global.end()) {
        if (json_global["brems"].is_string()) {
            std::string brems = json_global["brems"].get<std::string>();
            sec_def_global.utility_def.brems_def.parametrization
                = BremsstrahlungFactory::Get().GetEnumFromString(brems);
        } else {
            log_fatal("Invalid input for option 'brems'. Expected a string.");
        }
    } else {
        log_debug("The 'brems' option is not set. Use default %s",
            BremsstrahlungFactory::Get()
                .GetStringFromEnum(
                    sec_def_global.utility_def.brems_def.parametrization)
                .c_str());
    }

    if (json_global.find("ioniz") != json_global.end()) {
        if (json_global["ioniz"].is_string()) {
            std::string ioniz = json_global["ioniz"].get<std::string>();
            sec_def_global.utility_def.ioniz_def.parametrization
                = IonizationFactory::Get().GetEnumFromString(ioniz);
        } else {
            log_fatal("Invalid input for option 'ioniz'. Expected a string.");
        }
    } else {
        log_debug("The 'ioniz' option is not set. Use default %s",
            IonizationFactory::Get()
                .GetStringFromEnum(
                    sec_def_global.utility_def.ioniz_def.parametrization)
                .c_str());
    }

    if (json_global.find("epair") != json_global.end()) {
        if (json_global["epair"].is_string()) {
            std::string epair = json_global["epair"].get<std::string>();
            sec_def_global.utility_def.epair_def.parametrization
                = EpairProductionFactory::Get().GetEnumFromString(epair);
        } else {
            log_fatal("Invalid input for option 'epair'. Expected a string.");
        }
    } else {
        log_debug("The 'epair' option is not set. Use default %s",
            EpairProductionFactory::Get()
                .GetStringFromEnum(
                    sec_def_global.utility_def.epair_def.parametrization)
                .c_str());
    }

    if (json_global.find("photo") != json_global.end()) {
        if (json_global["photo"].is_string()) {
            std::string photo = json_global["photo"].get<std::string>();
            sec_def_global.utility_def.photo_def.parametrization
                = PhotonuclearFactory::Get().GetEnumFromString(photo);
        } else {
            log_fatal("Invalid input for option 'photo'. Expected a string.");
        }
    } else {
        log_debug("The 'photo' option is not set. Use default %s",
            PhotonuclearFactory::Get()
                .GetStringFromEnum(
                    sec_def_global.utility_def.photo_def.parametrization)
                .c_str());
    }

    if (json_global.find("annihilation") != json_global.end()) {
        if (json_global["annihilation"].is_string()) {
            std::string annihilation
                = json_global["annihilation"].get<std::string>();
            sec_def_global.utility_def.annihilation_def.parametrization
                = AnnihilationFactory::Get().GetEnumFromString(annihilation);
        } else {
            log_fatal(
                "Invalid input for option 'annihilation'. Expected a string.");
        }
    } else {
        log_debug("The 'annihilation' option is not set. Use default %s",
            AnnihilationFactory::Get()
                .GetStringFromEnum(
                    sec_def_global.utility_def.annihilation_def.parametrization)
                .c_str());
    }

    if (json_global.find("mupair") != json_global.end()) {
        if (json_global["mupair"].is_string()) {
            std::string mupair = json_global["mupair"].get<std::string>();
            sec_def_global.utility_def.mupair_def.parametrization
                = MupairProductionFactory::Get().GetEnumFromString(mupair);
        } else {
            log_fatal("Invalid input for option 'mupair'. Expected a string.");
        }
    } else {
        log_debug("The 'mupair' option is not set. Use default %s",
            MupairProductionFactory::Get()
                .GetStringFromEnum(
                    sec_def_global.utility_def.mupair_def.parametrization)
                .c_str());
    }

    if (json_global.find("weak") != json_global.end()) {
        if (json_global["weak"].is_string()) {
            std::string weak = json_global["weak"].get<std::string>();
            sec_def_global.utility_def.weak_def.parametrization
                = WeakInteractionFactory::Get().GetEnumFromString(weak);
        } else {
            log_fatal("Invalid input for option 'weak'. Expected a string.");
        }
    } else {
        log_debug("The 'weak' option is not set. Use default %s",
            WeakInteractionFactory::Get()
                .GetStringFromEnum(
                    sec_def_global.utility_def.weak_def.parametrization)
                .c_str());
    }

    if (json_global.find("compton") != json_global.end()) {
        if (json_global["compton"].is_string()) {
            std::string compton = json_global["compton"].get<std::string>();
            sec_def_global.utility_def.compton_def.parametrization
                = ComptonFactory::Get().GetEnumFromString(compton);
        } else {
            log_fatal("Invalid input for option 'compton'. Expected a string.");
        }
    } else {
        log_debug("The 'compton' option is not set. Use default %s",
            ComptonFactory::Get()
                .GetStringFromEnum(
                    sec_def_global.utility_def.compton_def.parametrization)
                .c_str());
    }

    if (json_global.find("photopair") != json_global.end()) {
        if (json_global["photopair"].is_string()) {
            std::string photopair = json_global["photopair"].get<std::string>();
            sec_def_global.utility_def.photopair_def.parametrization
                = PhotoPairFactory::Get().GetEnumFromString(photopair);
        } else {
            log_fatal(
                "Invalid input for option 'photopair'. Expected a string.");
        }
    } else {
        log_debug("The 'photopair' option is not set. Use default %s",
            PhotoPairFactory::Get()
                .GetStringFromEnum(
                    sec_def_global.utility_def.photopair_def.parametrization)
                .c_str());
    }

    if (json_global.find("photoangle") != json_global.end()) {
        if (json_global["photoangle"].is_string()) {
            std::string photoangle
                = json_global["photoangle"].get<std::string>();
            sec_def_global.utility_def.photopair_def.photoangle
                = PhotoPairFactory::Get().GetPhotoAngleEnumFromString(
                    photoangle);
        } else {
            log_fatal(
                "Invalid input for option 'photoangle'. Expected a string.");
        }
    } else {
        log_debug("The 'photoangle' option is not set. Use default %s",
            PhotoPairFactory::Get()
                .GetStringFromPhotoAngleEnum(
                    sec_def_global.utility_def.photopair_def.photoangle)
                .c_str());
    }

    if (json_global.find("photo_shadow") != json_global.end()) {
        if (json_global["photo_shadow"].is_string()) {
            std::string photo_shadow
                = json_global["photo_shadow"].get<std::string>();
            sec_def_global.utility_def.photo_def.shadow
                = PhotonuclearFactory::Get().GetShadowEnumFromString(
                    photo_shadow);
        } else {
            log_fatal(
                "Invalid input for option 'photo_shadow'. Expected a string.");
        }
    } else {
        log_debug("The 'photo_shadow' option is not set. Use default %s",
            PhotonuclearFactory::Get()
                .GetStringFromShadowEnum(
                    sec_def_global.utility_def.photo_def.shadow)
                .c_str());
    }

    if (json_global.find("photo_hard_component") != json_global.end()) {
        if (json_global["photo_hard_component"].is_boolean()) {
            sec_def_global.utility_def.photo_def.hard_component
                = json_global["photo_hard_component"].get<bool>();
        } else {
            log_fatal("Invalid input for option 'photo_hard_component'. "
                      "Expected a bool.");
        }
    } else {
        log_debug(
            "The 'photo_hard_component' option is not set. Use default %s",
            sec_def_global.utility_def.photo_def.hard_component ? "true"
                                                                : "false");
    }

    return sec_def_global;
}

std::string Propagator::ParseCutSettings(const std::string& json_object_str,
    const std::string& json_key, double default_ecut, double default_vcut,
    bool default_contrand)
{
    nlohmann::json json_object = nlohmann::json::parse(json_object_str);
    nlohmann::json output_object;
    std::string ecut_str = "e_cut";
    output_object[ecut_str] = default_ecut;
    std::string vcut_str = "v_cut";
    output_object[vcut_str] = default_vcut;
    std::string cont_str = "cont_rand";
    output_object[cont_str] = default_contrand;

    if (json_object.find(json_key) != json_object.end()) {
        if (json_object[json_key].is_object()) {
            json_object = json_object[json_key];
        } else {
            log_fatal("This cut setting is not a json object.");
        }
    } else {
        log_debug("The '%s' option is not set. Use default", json_key.c_str());
    }

    // get ecut
    if (json_object.find(ecut_str) != json_object.end()) {
        if (json_object[ecut_str].is_number()) {
            output_object[ecut_str] = json_object[ecut_str].get<double>();
        } else {
            log_fatal("Invalid input for option '%s'. Expected a number.",
                ecut_str.c_str());
        }
    } else {
        log_debug("The '%s' option is not set. Use default %f.",
            ecut_str.c_str(), default_ecut);
    }

    // get vcut
    if (json_object.find(vcut_str) != json_object.end()) {
        if (json_object[vcut_str].is_number()) {
            output_object[vcut_str] = json_object[vcut_str].get<double>();
        } else {
            log_fatal("Invalid input for option '%s'. Expected a number.",
                vcut_str.c_str());
        }
    } else {
        log_debug("The '%s' option is not set. Use default %f.",
            vcut_str.c_str(), default_vcut);
    }

    // get continuous randomization
    if (json_object.find(cont_str) != json_object.end()) {
        if (json_object[cont_str].is_boolean()) {
            output_object[cont_str] = json_object[cont_str].get<bool>();
        } else {
            log_fatal("Invalid input for option '%s'. Expected a number.",
                cont_str.c_str());
        }
    } else {
        log_debug("The '%s' option is not set. Use default %s.",
            cont_str.c_str(), default_contrand ? "true" : "false");
    }

    std::string output_object_str = output_object.dump();
    return output_object_str;
}
