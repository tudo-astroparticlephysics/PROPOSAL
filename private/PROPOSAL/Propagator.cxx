/*
 * Propagator.cxx
 *
 *  Created on: 23.04.2013
 *      Author: koehne
 */

// #include <cmath>

#include <fstream>

#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/geometry/Box.h"
#include "PROPOSAL/geometry/Cylinder.h"
#include "PROPOSAL/geometry/GeometryFactory.h"
#include "PROPOSAL/geometry/Sphere.h"

#include "PROPOSAL/medium/MediumFactory.h"

#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/json.hpp"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Global defaults
// ------------------------------------------------------------------------- //

const int Propagator::global_seed_            = 0;
const double Propagator::global_ecut_inside_  = 500;
const double Propagator::global_ecut_infront_ = -1;
const double Propagator::global_ecut_behind_  = -1;
const double Propagator::global_vcut_inside_  = -1;
const double Propagator::global_vcut_infront_ = 0.001;
const double Propagator::global_vcut_behind_  = -1;
const double Propagator::global_cont_inside_  = false;
const double Propagator::global_cont_infront_ = true;
const double Propagator::global_cont_behind_  = false;
const bool Propagator::do_interpolation_      = true;
const bool Propagator::uniform_               = true;

// ------------------------------------------------------------------------- //
// Constructors & destructor
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
Propagator::Propagator(const std::vector<Sector*>& sectors, const Geometry& geometry) try
    : current_sector_(NULL),
      particle_(sectors.at(0)->GetParticle()),
      detector_(geometry.clone())
{
    // --------------------------------------------------------------------- //
    // Check if all ParticleDefs are the same
    // --------------------------------------------------------------------- //

    for (auto sector: sectors)
    {
        if (sector->GetParticle().GetParticleDef() != particle_.GetParticleDef())
        {
            log_fatal("The particle definitions of the sectors must be identical for proper propagation!");
        } else
        {
            sectors_.push_back(new Sector(*sector));
        }
    }

    current_sector_ = sectors_.at(0);
} catch (const std::out_of_range& ex)
{
    log_fatal("No Sectors are provided for the Propagator!");
}

// ------------------------------------------------------------------------- //
Propagator::Propagator(const ParticleDef& particle_def,
                       const std::vector<Sector::Definition>& sector_defs,
                       const Geometry& geometry)
    : particle_(particle_def)
    , detector_(geometry.clone())
{
    for (auto def: sector_defs)
    {
        sectors_.push_back(new Sector(particle_, def));
    }

    try
    {
        current_sector_ = sectors_.at(0);
    } catch (const std::out_of_range& ex)
    {
        log_fatal("No Sectors are provided for the Propagator!");
    }
}

// ------------------------------------------------------------------------- //
Propagator::Propagator(const ParticleDef& particle_def,
                       const std::vector<Sector::Definition>& sector_defs,
                       const Geometry& geometry,
                       const InterpolationDef& interpolation_def)
    : particle_(particle_def)
    , detector_(geometry.clone())
{
    for (auto def: sector_defs)
    {
        sectors_.push_back(new Sector(particle_, def, interpolation_def));
    }

    try
    {
        current_sector_ = sectors_.at(0);
    } catch (const std::out_of_range& ex)
    {
        log_fatal("No Sectors are provided for the Propagator!");
    }
}

// ------------------------------------------------------------------------- //
Propagator::Propagator(const Propagator& propagator)
    : sectors_(propagator.sectors_.size(), NULL)
    , current_sector_(NULL)
    , particle_(propagator.particle_)
    , detector_(propagator.detector_->clone())
{
    for (unsigned int i = 0; i < propagator.sectors_.size(); ++i)
    {
        sectors_[i] = new Sector(particle_, *propagator.sectors_[i]);

        if (propagator.sectors_[i] == propagator.current_sector_)
        {
            current_sector_ = sectors_[i];
        }
    }
}

// ------------------------------------------------------------------------- //
Propagator::Propagator(const ParticleDef& particle_def, const std::string& config_file)
    : current_sector_(NULL)
    , particle_(particle_def)
    , detector_(NULL)
{
    int global_seed  = global_seed_;
    bool do_interpolation = do_interpolation_;
    bool uniform = uniform_;

    InterpolationDef interpolation_def;
    Sector::Definition sec_def_global;

    // Create the json parser
    nlohmann::json json_config;
    try
    {
        std::string expanded_config_file_path = Helper::ResolvePath(config_file);
        std::ifstream infilestream(expanded_config_file_path);
        infilestream >> json_config;
    }
    catch (const nlohmann::json::parse_error& e)
    {
        log_fatal("Unable parse \"%s\" as json file", config_file.c_str());
    }

    // set global settings
    if (json_config.find("global") == json_config.end())
    {
        log_fatal("No given global settings. Use default");
    }
    if (!json_config["global"].is_object())
    {
        log_fatal("The given global option is not an object.");
    }
    nlohmann::json json_global = json_config["global"];
    std::string json_global_str = json_global.dump();

    // set the seed of the random number generator
    if (json_global.find("seed") != json_global.end())
    {
        if (json_global["seed"].is_number())
        {
            global_seed = json_global["seed"].get<int>();
            // this can throw the error (nlohmann::detail::type_error& e) if no numbercheck
        }
        else
        {
            log_fatal("The given seed option is not a number.");
        }
    }

    RandomGenerator::Get().SetSeed(global_seed);
    log_info("Seed of the default random generator set to %i", global_seed);

    // Read in global cut and continous randomization options
    std::string cut_object_str;
    nlohmann::json cuts_infront_object;
    cut_object_str = ParseCutSettings(json_global_str,
                                      "cuts_infront",
                                      global_ecut_infront_,
                                      global_vcut_infront_,
                                      global_cont_infront_);
    cuts_infront_object = nlohmann::json::parse(cut_object_str);

    nlohmann::json cuts_inside_object;
    cut_object_str = ParseCutSettings(json_global_str,
                                      "cuts_inside",
                                      global_ecut_inside_,
                                      global_vcut_inside_,
                                      global_cont_inside_);
    cuts_inside_object = nlohmann::json::parse(cut_object_str);

    nlohmann::json cuts_behind_object;
    cut_object_str = ParseCutSettings(json_global_str,
                                      "cuts_behind",
                                      global_ecut_behind_,
                                      global_vcut_behind_,
                                      global_cont_behind_);
    cuts_behind_object = nlohmann::json::parse(cut_object_str);

    // Read in interpolation options
    if (json_global.find("interpolation") != json_global.end())
    {
        if (json_global["interpolation"].is_object())
        {
            if (json_global["interpolation"].find("do_interpolation") != json_global["interpolation"].end())
            {
                if (json_global["interpolation"]["do_interpolation"].is_boolean())
                {
                    do_interpolation = json_global["interpolation"]["do_interpolation"];
                }
                else
                {
                    log_fatal("The given do_interpolation option is not a bool.");
                }
            }
            else
            {
                log_debug("The do_interpolation option is not set. Default is true");
            }
            if (do_interpolation)
            {
                std::string json_str = json_global["interpolation"].dump();
                interpolation_def = CreateInterpolationDef(json_str);
            }
            else
            {
                log_info("Integration instead of interpolation is chosen. The propagation is now extreamly slow.");
            }
        }
        else
        {
            log_fatal("The given global interpolation option is not a json object.");
        }
    }
    else
    {
        log_debug("No Interpolation-Option set. Use default InterpolationDef");
    }

    // Read in global sector definition
    sec_def_global = CreateSectorDefinition(json_global_str);

    // read in the uniform flag
    if (json_global.find("uniform") != json_global.end())
    {
        if (json_global["uniform"].is_boolean())
        {
            uniform = json_global["uniform"].get<bool>();
        }
        else
        {
            log_fatal("The given uniform phasespace sampling option for decays is not a bool.");
        }
    }
    else
    {
        log_debug("No uniform phasespace sampling option for decays set. Use default (true)");
    }
    particle_.GetDecayTable().SetUniformSampling(uniform);


    // Parse detector geometry
    if (json_config.find("detector") != json_config.end())
    {
        std::string json_str = json_config["detector"].dump();
        detector_ = ParseGeometryConifg(json_str);
    }
    else
    {
        log_fatal("You need to specify a detector geometry.");
    }

    // Read in all sector definitions
    nlohmann::json json_sectors;

    if (json_config.find("sectors") != json_config.end())
    {
        if (json_config["sectors"].is_array())
        {
            if (1 > json_config["sectors"].size())
            {
                log_fatal("You need to specify at least one sector.");
            }
        }
        else
        {
            log_fatal("The given sectors option is not an array of json objects.");
        }
    }
    else
    {
        log_fatal("No given sectors option.");
    }


    for (size_t idx=0; idx < json_config["sectors"].size(); idx++)
    {
        json_sectors = json_config["sectors"][idx];
        if (!json_sectors.is_object())
        {
            log_fatal("This given sector is not a json object.");
        }
        std::string json_sector_str = json_sectors.dump();

        // Create medium
        MediumFactory::Definition medium_def;
        if (json_sectors.find("density_correction") != json_sectors.end())
        {
            if (json_sectors["density_correction"].is_number())
            {
                medium_def.density_correction = json_sectors["density_correction"].get<double>();
            }
            else
            {
                log_fatal("The given density correction option is not a number.");
            }
        }
        else
        {
            log_debug("No given density correction. Use default 1.0");
        }

        std::string medium_name = MediumFactory::Get().GetStringFromEnum(medium_def.type);
        if (json_sectors.find("medium") != json_sectors.end())
        {
            if (json_sectors["medium"].is_string())
            {
                medium_name = json_sectors["medium"].get<std::string>();
            }
            else
            {
                log_fatal("The given medium option is not a string.");
            }
        }
        else
        {
            log_debug("No given medium. Use default (Water)");
        }
        Medium* med = MediumFactory::Get().CreateMedium(medium_name, medium_def.density_correction);

        // Create Geometry
        Geometry* geometry = NULL;

        if (json_sectors.find("geometry") != json_sectors.end())
        {
            std::string json_str = json_sectors["geometry"].dump();
            geometry = ParseGeometryConifg(json_str);
        }
        else
        {
            log_fatal("You need to specify a geometry for each sector");
        }

        double hierarchy = geometry->GetHierarchy();
        if (json_sectors.find("hierarchy") != json_sectors.end())
        {
            if (json_sectors["hierarchy"].is_number())
            {
                hierarchy = json_sectors["hierarchy"].get<int>();
            }
            else
            {
                log_fatal("The given hierarchy option is not a number.");
            }
        }
        else
        {
            log_debug("No given hierarchy. Use default (0)");
        }
        geometry->SetHierarchy(hierarchy);

        // Use global options in case they will not be overriden
        Sector::Definition sec_def_infront = sec_def_global;
        sec_def_infront.location           = Sector::ParticleLocation::InfrontDetector;
        sec_def_infront.SetMedium(*med);
        sec_def_infront.SetGeometry(*geometry);

        Sector::Definition sec_def_inside = sec_def_global;
        sec_def_inside.location           = Sector::ParticleLocation::InsideDetector;
        sec_def_inside.SetMedium(*med);
        sec_def_inside.SetGeometry(*geometry);

        Sector::Definition sec_def_behind = sec_def_global;
        sec_def_behind.location           = Sector::ParticleLocation::BehindDetector;
        sec_def_behind.SetMedium(*med);
        sec_def_behind.SetGeometry(*geometry);

        nlohmann::json json_cutsettings;

        // cut settings infront
        cut_object_str = ParseCutSettings(json_sector_str,
                                          "cuts_infront",
                                          cuts_infront_object["e_cut"].get<double>(),
                                          cuts_infront_object["v_cut"].get<double>(),
                                          cuts_infront_object["cont_rand"].get<bool>());

        json_cutsettings = nlohmann::json::parse(cut_object_str);
        sec_def_infront.cut_settings.SetEcut(json_cutsettings["e_cut"].get<double>());
        sec_def_infront.cut_settings.SetVcut(json_cutsettings["v_cut"].get<double>());
        sec_def_infront.do_continuous_randomization = json_cutsettings["cont_rand"].get<bool>();

        // cut settings inside
        cut_object_str = ParseCutSettings(json_sector_str,
                                          "cuts_inside",
                                          cuts_inside_object["e_cut"].get<double>(),
                                          cuts_inside_object["v_cut"].get<double>(),
                                          cuts_inside_object["cont_rand"].get<bool>());

        json_cutsettings = nlohmann::json::parse(cut_object_str);
        sec_def_inside.cut_settings.SetEcut(json_cutsettings["e_cut"].get<double>());
        sec_def_inside.cut_settings.SetVcut(json_cutsettings["v_cut"].get<double>());
        sec_def_inside.do_continuous_randomization = json_cutsettings["cont_rand"].get<bool>();

        // cut settings behind
        cut_object_str = ParseCutSettings(json_sector_str,
                                          "cuts_behind",
                                          cuts_behind_object["e_cut"].get<double>(),
                                          cuts_behind_object["v_cut"].get<double>(),
                                          cuts_behind_object["cont_rand"].get<bool>());

        json_cutsettings = nlohmann::json::parse(cut_object_str);
        sec_def_behind.cut_settings.SetEcut(json_cutsettings["e_cut"].get<double>());
        sec_def_behind.cut_settings.SetVcut(json_cutsettings["v_cut"].get<double>());
        sec_def_behind.do_continuous_randomization = json_cutsettings["cont_rand"].get<bool>();

        if (do_interpolation)
        {
            sectors_.push_back(new Sector(particle_, sec_def_infront, interpolation_def));
            sectors_.push_back(new Sector(particle_, sec_def_inside, interpolation_def));
            sectors_.push_back(new Sector(particle_, sec_def_behind, interpolation_def));
        } else
        {
            sectors_.push_back(new Sector(particle_, sec_def_infront));
            sectors_.push_back(new Sector(particle_, sec_def_inside));
            sectors_.push_back(new Sector(particle_, sec_def_behind));
        }

        delete geometry;
        delete med;
    }
}

Propagator::~Propagator()
{
    for (auto sector: sectors_)
    {
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
    if (*detector_ != *propagator.detector_)
    {
        return false;
    }
    if (sectors_.size() != propagator.sectors_.size())
    {
        return false;
    }
    else
    {
        for (unsigned int i = 0; i < sectors_.size(); ++i)
        {
            if (*sectors_[i] != *propagator.sectors_[i])
            {
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
std::vector<DynamicData*> Propagator::Propagate(double MaxDistance_cm)
{
    Output::getInstance().ClearSecondaryVector();

    Output::getInstance().GetSecondarys().reserve(1000);

#if ROOT_SUPPORT
    Output::getInstance().StorePrimaryInTree(&particle_);
#endif

    if (Output::store_in_ASCII_file_)
        Output::getInstance().StorePrimaryInASCII(&particle_);

    double distance = 0;
    double result   = 0;

    // These two variables are needed to calculate the energy loss inside the detector
    // energy_at_entry_point is initialized with the current energy because this is a
    // reasonable value for particle_ which starts inside the detector

    double energy_at_entry_point = particle_.GetEnergy();
    double energy_at_exit_point  = 0;

    Vector3D particle_position  = particle_.GetPosition();
    Vector3D particle_direction = particle_.GetDirection();

    bool starts_in_detector = detector_->IsInside(particle_position, particle_direction);
    bool is_in_detector     = false;
    bool was_in_detector    = false;

    while (1)
    {
        particle_position  = particle_.GetPosition();
        particle_direction = particle_.GetDirection();

        ChooseCurrentCollection(particle_position, particle_direction);

        if (current_sector_ == NULL)
        {
            log_info("particle_ reached the border");
            break;
        }

        // Check if have to propagate the particle_ through the whole collection
        // or only to the collection border
        distance = CalculateEffectiveDistance(particle_position, particle_direction);

        is_in_detector = detector_->IsInside(particle_position, particle_direction);
        // entry point of the detector
        if (!starts_in_detector && !was_in_detector && is_in_detector)
        {
            particle_.SetEntryPoint(particle_position);
            particle_.SetEntryEnergy(particle_.GetEnergy());
            particle_.SetEntryTime(particle_.GetTime());

            energy_at_entry_point = particle_.GetEnergy();

            was_in_detector = true;
        }
        // exit point of the detector
        else if (was_in_detector && !is_in_detector)
        {
            particle_.SetExitPoint(particle_position);
            particle_.SetExitEnergy(particle_.GetEnergy());
            particle_.SetExitTime(particle_.GetTime());

            energy_at_exit_point = particle_.GetEnergy();
            // we don't want to run in this case a second time so we set was_in_detector to false
            was_in_detector = false;

        }
        // if particle_ starts inside the detector we only ant to fill the exit point
        else if (starts_in_detector && !is_in_detector)
        {
            particle_.SetExitPoint(particle_position);
            particle_.SetExitEnergy(particle_.GetEnergy());
            particle_.SetExitTime(particle_.GetTime());

            energy_at_exit_point = particle_.GetEnergy();
            // we don't want to run in this case a second time so we set starts_in_detector to false
            starts_in_detector = false;
        }
        if (MaxDistance_cm <= particle_.GetPropagatedDistance() + distance)
        {
            distance = MaxDistance_cm - particle_.GetPropagatedDistance();
        }

        result = current_sector_->Propagate(distance);

        if (result <= 0 || MaxDistance_cm <= particle_.GetPropagatedDistance())
            break;
    }

    particle_.SetElost(energy_at_entry_point - energy_at_exit_point);

#if ROOT_SUPPORT
    Output::getInstance().StorePropagatedPrimaryInTree(&particle_);
#endif
    if (Output::store_in_ASCII_file_)
        Output::getInstance().StorePropagatedPrimaryInASCII(&particle_);

    return Output::getInstance().GetSecondarys();
}

// ------------------------------------------------------------------------- //
void Propagator::ChooseCurrentCollection(const Vector3D& particle_position, const Vector3D& particle_direction)
{
    std::vector<int> crossed_collections;
    crossed_collections.resize(0);

    for (unsigned int i = 0; i < sectors_.size(); i++)
    {
        if (detector_->IsInfront(particle_position, particle_direction))
        {
            if (sectors_[i]->GetLocation() != 0)
            {
                continue;
            } else
            {
                if (sectors_[i]->GetGeometry()->IsInside(particle_position, particle_direction))
                {
                    current_sector_ = sectors_[i];
                    crossed_collections.push_back(i);
                } else
                {
                }
            }
        }

        else if (detector_->IsInside(particle_position, particle_direction))
        {
            if (sectors_[i]->GetLocation() != 1)
            {
                continue;
            } else
            {
                if (sectors_[i]->GetGeometry()->IsInside(particle_position, particle_direction))
                {
                    current_sector_ = sectors_[i];
                    crossed_collections.push_back(i);
                } else
                {
                }
            }

        }

        else if (detector_->IsBehind(particle_position, particle_direction))
        {
            if (sectors_[i]->GetLocation() != 2)
            {
                continue;
            } else
            {
                if (sectors_[i]->GetGeometry()->IsInside(particle_position, particle_direction))
                {
                    current_sector_ = sectors_[i];
                    crossed_collections.push_back(i);
                }
                // The particle reached the border of all specified collections
                else
                {
                }
            }
        }
    }

    // No process collection was found
    if (crossed_collections.size() == 0)
    {
        current_sector_ = NULL;
        log_fatal("No Cross Section was found!!!");
    }

    // Choose current collection when multiple collections are crossed!
    //
    // Choose by hierarchy of Geometry
    // If same hierarchys are available the denser one is choosen
    // If hierarchy and density are the same then the first found is taken.
    //

    for (std::vector<int>::iterator iter = crossed_collections.begin(); iter != crossed_collections.end(); ++iter)
    {
        // Current Hierachy is bigger -> Nothing to do!
        //
        if (current_sector_->GetGeometry()->GetHierarchy() > sectors_[*iter]->GetGeometry()->GetHierarchy())
        {
            continue;
        }
        // Current Hierachy is equal -> Look at the density!
        //
        else if (current_sector_->GetGeometry()->GetHierarchy() == sectors_[*iter]->GetGeometry()->GetHierarchy())
        {
            // Current Density is bigger or same -> Nothing to do!
            //

            if (current_sector_->GetMedium()->GetMassDensity() >= sectors_[*iter]->GetMedium()->GetMassDensity())
            {
                continue;
            }

            // Current Density is smaller -> Set the new collection!
            //
            else
            {
                current_sector_ = sectors_[*iter];
            }

        }

        // Current Hierachy is smaller -> Set the new collection!
        //
        else
        {
            current_sector_ = sectors_[*iter];
        }
    }
}

// ------------------------------------------------------------------------- //
double Propagator::CalculateEffectiveDistance(const Vector3D& particle_position, const Vector3D& particle_direction)
{
    double distance_to_collection_border = 0;
    double distance_to_detector          = 0;
    double distance_to_closest_approach  = 0;

    distance_to_collection_border =
        current_sector_->GetGeometry()->DistanceToBorder(particle_position, particle_direction).first;
    double tmp_distance_to_border;

    for (std::vector<Sector*>::iterator iter = sectors_.begin(); iter != sectors_.end(); ++iter)
    {
        if (detector_->IsInfront(particle_position, particle_direction))
        {
            if ((*iter)->GetLocation() != 0)
                continue;
            else
            {
                if ((*iter)->GetGeometry()->GetHierarchy() >= current_sector_->GetGeometry()->GetHierarchy())
                {
                    tmp_distance_to_border =
                        (*iter)->GetGeometry()->DistanceToBorder(particle_position, particle_direction).first;
                    if (tmp_distance_to_border <= 0)
                        continue;
                    distance_to_collection_border = std::min(tmp_distance_to_border, distance_to_collection_border);
                }
            }
        }

        else if (detector_->IsInside(particle_position, particle_direction))
        {
            if ((*iter)->GetLocation() != 1)
                continue;
            else
            {
                tmp_distance_to_border =
                    (*iter)->GetGeometry()->DistanceToBorder(particle_position, particle_direction).first;
                if (tmp_distance_to_border <= 0)
                    continue;
                distance_to_collection_border = std::min(tmp_distance_to_border, distance_to_collection_border);
            }

        }

        else if (detector_->IsBehind(particle_position, particle_direction))
        {
            if ((*iter)->GetLocation() != 2)
                continue;
            else
            {
                if ((*iter)->GetGeometry()->GetHierarchy() >= current_sector_->GetGeometry()->GetHierarchy())
                {
                    tmp_distance_to_border =
                        (*iter)->GetGeometry()->DistanceToBorder(particle_position, particle_direction).first;
                    if (tmp_distance_to_border <= 0)
                        continue;
                    distance_to_collection_border = std::min(tmp_distance_to_border, distance_to_collection_border);
                }
                // The particle_ reached the border of all specified collections
                else
                {
                }
            }
        }
    }

    distance_to_detector = detector_->DistanceToBorder(particle_position, particle_direction).first;

    distance_to_closest_approach = detector_->DistanceToClosestApproach(particle_position, particle_direction);

    if (std::abs(distance_to_closest_approach) < GEOMETRY_PRECISION)
    {
        particle_.SetClosestApproachPoint(particle_position);
        particle_.SetClosestApproachEnergy(particle_.GetEnergy());
        particle_.SetClosestApproachTime(particle_.GetTime());

        distance_to_closest_approach = 0;
    }

    if (distance_to_detector > 0)
    {
        if (distance_to_closest_approach > 0)
        {
            return std::min({distance_to_detector, distance_to_collection_border, distance_to_closest_approach});
        } else
        {
            return std::min(distance_to_detector, distance_to_collection_border);
        }
    } else
    {
        if (distance_to_closest_approach > 0)
        {
            return std::min(distance_to_closest_approach, distance_to_collection_border);
        } else
        {
            return distance_to_collection_border;
        }
    }
}

// ------------------------------------------------------------------------- //
// Private member functions
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
Geometry* Propagator::ParseGeometryConifg(const std::string& json_object_str)
{
    nlohmann::json json_object = nlohmann::json::parse(json_object_str);
    if (!json_object.is_object())
    {
        log_fatal("This geometry option should be a json object.");
    }

    std::string origin_str       = "origin";
    std::string outer_radius_str = "outer_radius";
    std::string inner_radius_str = "inner_radius";
    std::string lenght_str       = "lenght";
    std::string width_str        = "width";
    std::string height_str       = "height";

    std::string warning_str = "Geometry %s needs to specify \"%s\" in the config file!";

    double cm_to_meter = 100.0;

    // --------------------------------------------------------------------- //
    // Get geometry from default constructor
    // --------------------------------------------------------------------- //

    Geometry* geometry = NULL;

    if (json_object.find("shape") != json_object.end())
    {
        if (json_object["shape"].is_string())
        {
            std::string name = json_object["shape"].get<std::string>();
            geometry = GeometryFactory::Get().CreateGeometry(name);
        }
        else
        {
            log_fatal("The given shape option is not a string.");
        }
    }
    else
    {
        log_fatal(warning_str.c_str(), "", "shape");
    }

    // --------------------------------------------------------------------- //
    // Get the position vector from the property tree
    // --------------------------------------------------------------------- //

    double x = 0;
    double y = 0;
    double z = 0;

    if (json_object.find(origin_str) != json_object.end())
    {
        if (json_object[origin_str].is_array())
        {
            if (json_object[origin_str].size() == 3)
            {
                for (size_t idx=0; idx < json_object[origin_str].get<std::vector<double>>().size(); idx++)
                {
                    if (json_object[origin_str][idx].is_number())
                    {
                        double coord = json_object[origin_str][idx].get<double>() * cm_to_meter;

                        switch (idx)
                        {
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
                    }
                    else
                    {
                        log_fatal("The given origin option needs doubles in the array.");
                    }
                }
            }
            else
            {
                log_fatal("The given origin option is not an array of len 3.");
            }
        }
        else
        {
            log_fatal("The given origin option is not an array of 3 doubles.");
        }
    }
    else
    {
        log_fatal(warning_str.c_str(), "", "origin");
    }

    Vector3D vec(x, y, z);

    // --------------------------------------------------------------------- //
    // Check type of geometry and set members
    // --------------------------------------------------------------------- //

    if (PROPOSAL::Sphere* sphere = dynamic_cast<PROPOSAL::Sphere*>(geometry))
    {
        double radius = 0;
        double inner_radius = 0;

        if (json_object.find(outer_radius_str) != json_object.end())
        {
            if (json_object[outer_radius_str].is_number())
            {
                radius = json_object[outer_radius_str].get<double>() * cm_to_meter;
            }
            else
            {
                log_fatal("The given outer_radius option is not a number.");
            }
        }
        else
        {
            log_fatal(warning_str.c_str(), sphere->GetName().c_str(), outer_radius_str.c_str());
        }

        if (json_object.find(inner_radius_str) != json_object.end())
        {
            if (json_object[inner_radius_str].is_number())
            {
                inner_radius = json_object[inner_radius_str].get<double>() * cm_to_meter;
            }
            else
            {
                log_fatal("The given inner_radius option is not a number.");
            }
        }
        else
        {
            log_fatal(warning_str.c_str(), sphere->GetName().c_str(), inner_radius_str.c_str());
        }


        sphere->SetPosition(vec);
        sphere->SetRadius(radius);
        sphere->SetInnerRadius(inner_radius);

        return sphere;
    } else if (PROPOSAL::Box* box = dynamic_cast<PROPOSAL::Box*>(geometry))
    {
        double x = 0;
        double y = 0;
        double z = 0;

        if (json_object.find(lenght_str) != json_object.end())
        {
            if (json_object[lenght_str].is_number())
            {
                x = json_object[lenght_str].get<double>() * cm_to_meter;
            }
            else
            {
                log_fatal("The given length option is not a number.");
            }
        }
        else
        {
            log_fatal(warning_str.c_str(), sphere->GetName().c_str(), lenght_str.c_str());
        }

        if (json_object.find(width_str) != json_object.end())
        {
            if (json_object[width_str].is_number())
            {
                y = json_object[width_str].get<double>() * cm_to_meter;
            }
            else
            {
                log_fatal("The given width option is not a number.");
            }
        }
        else
        {
            log_fatal(warning_str.c_str(), sphere->GetName().c_str(), width_str.c_str());
        }

        if (json_object.find(height_str) != json_object.end())
        {
            if (json_object[height_str].is_number())
            {
                z = json_object[height_str].get<double>() * cm_to_meter;
            }
            else
            {
                log_fatal("The given height option is not a number.");
            }
        }
        else
        {
            log_fatal(warning_str.c_str(), sphere->GetName().c_str(), height_str.c_str());
        }


        box->SetPosition(vec);
        box->SetX(x);
        box->SetY(y);
        box->SetZ(z);

        return box;
    } else if (PROPOSAL::Cylinder* cylinder = dynamic_cast<PROPOSAL::Cylinder*>(geometry))
    {
        double radius = 0;
        double inner_radius = 0;
        double z = 0;

        if (json_object.find(outer_radius_str) != json_object.end())
        {
            if (json_object[outer_radius_str].is_number())
            {
                radius = json_object[outer_radius_str].get<double>() * cm_to_meter;
            }
            else
            {
                log_fatal("The given outer_radius option is not a number.");
            }
        }
        else
        {
            log_fatal(warning_str.c_str(), sphere->GetName().c_str(), outer_radius_str.c_str());
        }

        if (json_object.find(inner_radius_str) != json_object.end())
        {
            if (json_object[inner_radius_str].is_number())
            {
                inner_radius = json_object[inner_radius_str].get<double>() * cm_to_meter;
            }
            else
            {
                log_fatal("The given inner_radius option is not a number.");
            }
        }
        else
        {
            log_fatal(warning_str.c_str(), sphere->GetName().c_str(), inner_radius_str.c_str());
        }

        if (json_object.find(height_str) != json_object.end())
        {
            if (json_object[height_str].is_number())
            {
                z = json_object[height_str].get<double>() * cm_to_meter;
            }
            else
            {
                log_fatal("The given height option is not a number.");
            }
        }
        else
        {
            log_fatal(warning_str.c_str(), sphere->GetName().c_str(), height_str.c_str());
        }


        cylinder->SetPosition(vec);
        cylinder->SetRadius(radius);
        cylinder->SetInnerRadius(inner_radius);
        cylinder->SetZ(z);

        return cylinder;
    } else
    {
        log_fatal("Dynamic casts of Geometries failed. Should not end here!");
        return NULL;
    }
}

InterpolationDef Propagator::CreateInterpolationDef(const std::string& json_object_str)
{
    // Read in interpolation options
    InterpolationDef interpolation_def;
    nlohmann::json json_object = nlohmann::json::parse(json_object_str);

    if (json_object.find("do_binary_tables") != json_object.end())
    {
        if (json_object["do_binary_tables"].is_boolean())
        {
            interpolation_def.raw = json_object["do_binary_tables"];
        }
        else
        {
            log_fatal("The given do_binary_tables option is not a bool.");
        }
    }
    else
    {
        log_debug("The do_binary_tables option is not set. Use default (true");
    }

    // Parse to find path to interpolation tables
    if (json_object.find("path_to_tables") != json_object.end())
    {
        std::string table_path_str = "";
        if (json_object["path_to_tables"].is_string())
        {
            table_path_str = json_object["path_to_tables"];
            interpolation_def.path_to_tables = Helper::ResolvePath(table_path_str);
        }
        else if (json_object["path_to_tables"].is_array())
        {
            for (size_t idx=0; idx < json_object["path_to_tables"].get<std::vector<std::string>>().size(); idx++)
            {
                if (json_object["path_to_tables"][idx].is_string())
                {
                    table_path_str = Helper::ResolvePath(json_object["path_to_tables"][idx].get<std::string>());
                    if (table_path_str != "")
                        break;
                }
                else
                {
                    log_fatal("The given path_to_tables option does not consist of strings.");
                }
            }
        }
        else
        {
            log_fatal("The given path_to_tables option must be a string or a list of strings.");
        }

        if (table_path_str != "")
        {
            interpolation_def.path_to_tables = table_path_str;
            log_info("Path to interpolation tables set to: \"%s\"", table_path_str.c_str());
        }
        else
        {
            log_warn("No valid path to interpolation tables found. Save tables in memory!");
        }
    }
    else
    {
        log_debug("No path to tables set. Use default and save in memory");
    }

    return interpolation_def;
}

Sector::Definition Propagator::CreateSectorDefinition(const std::string& json_object_str)
{
    Sector::Definition sec_def_global;
    nlohmann::json json_global = nlohmann::json::parse(json_object_str);

    if (json_global.find("brems_multiplier") != json_global.end())
    {
        if (json_global["brems_multiplier"].is_number())
        {
            sec_def_global.utility_def.brems_def.multiplier = json_global["brems_multiplier"].get<double>();
        }
        else
        {
            log_fatal("The given brems_multiplier option is not a double.");
        }
    }
    else
    {
        log_debug("No given brems_multiplier option given. Use default (%f)", sec_def_global.utility_def.brems_def.multiplier);
    }

    if (json_global.find("photo_multiplier") != json_global.end())
    {
        if (json_global["photo_multiplier"].is_number())
        {
            sec_def_global.utility_def.photo_def.multiplier = json_global["photo_multiplier"].get<double>();
        }
        else
        {
            log_fatal("The given photo_multiplier option is not a double.");    
        }
        
    }
    else
    {
        log_debug("No given photo_multiplier option given. Use default (%f)", sec_def_global.utility_def.photo_def.multiplier);
    }

    if (json_global.find("ioniz_multiplier") != json_global.end())
    {
        if (json_global["ioniz_multiplier"].is_number())
        {
            sec_def_global.utility_def.ioniz_def.multiplier = json_global["ioniz_multiplier"].get<double>();
        }
        else
        {
            log_fatal("The given ioniz_multiplier option is not a double.");
        }
    }
    else
    {
        log_debug("No given ioniz_multiplier option given. Use default (%f)", sec_def_global.utility_def.ioniz_def.multiplier);
    }

    if (json_global.find("epair_multiplier") != json_global.end())
    {
        if (json_global["epair_multiplier"].is_number())
        {
            sec_def_global.utility_def.epair_def.multiplier = json_global["epair_multiplier"].get<double>();
        }
        else
        {
            log_fatal("The given epair_multiplier option is not a double.");
        }
    }
    else
    {
        log_debug("No given epair_multiplier option given. Use default (%f)", sec_def_global.utility_def.epair_def.multiplier);
    }

    if (json_global.find("lpm") != json_global.end())
    {
        if (json_global["lpm"].is_boolean())
        {
            sec_def_global.utility_def.epair_def.lpm_effect = json_global["lpm"].get<bool>();
            sec_def_global.utility_def.brems_def.lpm_effect = json_global["lpm"].get<bool>();
        }
        else
        {
            log_fatal("The given lpm option is not a bool.");
        }
    }
    else
    {
        log_debug("No given lpm option given. Use default (true)");
    }

    if (json_global.find("exact_time") != json_global.end())
    {
        if (json_global["exact_time"].is_boolean())
        {
            sec_def_global.do_exact_time_calculation = json_global["exact_time"].get<bool>();
        }
        else
        {
            log_fatal("The given exact_time option is not a bool.");
        }
    }
    else
    {
        log_debug("No given exact_time option given. Use default (true)");
    }

    if (json_global.find("continous_loss_output") != json_global.end())
    {
        if (json_global["continous_loss_output"].is_boolean())
        {
            sec_def_global.do_continuous_energy_loss_output = json_global["continous_loss_output"].get<bool>();
        }
        else
        {
            log_fatal("The given continous_loss_output option is not a bool.");
        }
    }
    else
    {
        log_debug("No given continous_loss_output option given. Use default true");
    }

    if (json_global.find("stopping_decay") != json_global.end())
    {
        if (json_global["stopping_decay"].is_boolean())
        {
            sec_def_global.stopping_decay = json_global["stopping_decay"].get<bool>();
        }
        else
        {
            log_fatal("The given stopping_decay option is not a bool.");
        }
    }
    else
    {
        log_debug("No given stopping_decay option given. Use default true");
    }

    if (json_global.find("stochastic_loss_weighting") != json_global.end())
    {
        if (json_global["stochastic_loss_weighting"].is_number())
        {
            sec_def_global.stochastic_loss_weighting = json_global["stochastic_loss_weighting"].get<double>();
            sec_def_global.do_stochastic_loss_weighting = 1e-10 < std::abs(sec_def_global.stochastic_loss_weighting);
        }
        else
        {
            log_fatal("The given stochastic_loss_weighting option is not a number.");
        }
    }
    else
    {
        log_debug("The stochastic_loss_weighting option is not set. Use default (%f)", sec_def_global.stochastic_loss_weighting);
    }


    if (json_global.find("scattering") != json_global.end())
    {
        if (json_global["scattering"].is_string())
        {
            std::string scattering = json_global["scattering"].get<std::string>();
            sec_def_global.scattering_model = ScatteringFactory::Get().GetEnumFromString(scattering);
        }
        else
        {
            log_fatal("The given scattering option is not a string.");
        }
    }
    else
    {
        log_debug("The scattering option is not set. Use default (%s)",
                  ScatteringFactory::Get().GetStringFromEnum(sec_def_global.scattering_model).c_str());
    }

    if (json_global.find("brems") != json_global.end())
    {
        if (json_global["brems"].is_string())
        {
            std::string brems = json_global["brems"].get<std::string>();
            sec_def_global.utility_def.brems_def.parametrization = BremsstrahlungFactory::Get().GetEnumFromString(brems);
        }
        else
        {
            log_fatal("The given brems option is not a string.");
        }
    }
    else
    {
        log_debug("The brems option is not set. Use default %s",
                  BremsstrahlungFactory::Get().GetStringFromEnum(sec_def_global.utility_def.brems_def.parametrization).c_str());
    }

    if (json_global.find("epair") != json_global.end())
    {
        if (json_global["epair"].is_string())
        {
            std::string epair = json_global["epair"].get<std::string>();
            sec_def_global.utility_def.epair_def.parametrization = EpairProductionFactory::Get().GetEnumFromString(epair);
        }
        else
        {
            log_fatal("The given epair option is not a string.");
        }
    }
    else
    {
        log_debug("The epair option is not set. Use default %s",
                  EpairProductionFactory::Get().GetStringFromEnum(sec_def_global.utility_def.epair_def.parametrization).c_str());
    }
    
    if (json_global.find("photo") != json_global.end())
    {
        if (json_global["photo"].is_string())
        {
            std::string photo = json_global["photo"].get<std::string>();
            sec_def_global.utility_def.photo_def.parametrization = PhotonuclearFactory::Get().GetEnumFromString(photo);
        }
        else
        {
            log_fatal("The given photo option is not a string.");
        }
    }
    else
    {
        log_debug("The photo option is not set. Use default %s",
                  PhotonuclearFactory::Get().GetStringFromEnum(sec_def_global.utility_def.photo_def.parametrization).c_str());
    }

    if (json_global.find("photo_shadow") != json_global.end())
    {
        if (json_global["photo_shadow"].is_string())
        {
            std::string photo_shadow = json_global["photo_shadow"].get<std::string>();
            sec_def_global.utility_def.photo_def.shadow = PhotonuclearFactory::Get().GetShadowEnumFromString(photo_shadow);
        }
        else
        {
            log_fatal("The given photo_shadow option is not a string.");
        }
    }
    else
    {
        log_debug("The photo_shadow option is not set. Use default %s",
                  PhotonuclearFactory::Get().GetStringFromShadowEnum(sec_def_global.utility_def.photo_def.shadow).c_str());
    }

    if (json_global.find("photo_hard_component") != json_global.end())
    {
        if (json_global["photo_hard_component"].is_boolean())
        {
            sec_def_global.utility_def.photo_def.hard_component = json_global["photo_hard_component"].get<bool>();
        }
        else
        {
            log_fatal("The given photo_hard_component option is not a bool.");
        }
    }
    else
    {
        log_debug("The photo_hard_component option is not set. Use default %s",
                  sec_def_global.utility_def.photo_def.hard_component ? "true" : "false");
    }

    return sec_def_global;
}


std::string Propagator::ParseCutSettings(const std::string& json_object_str,
                                         const std::string& json_key,
                                         double default_ecut,
                                         double default_vcut,
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

    if (json_object.find(json_key) != json_object.end())
    {
        if (json_object[json_key].is_object())
        {
            json_object = json_object[json_key];
        }
        else
        {
            log_fatal("This cut setting is not a json object.");
        }
    }
    else
    {
        log_debug("No given %s . Use default", json_key.c_str());
    }

    // get ecut
    if (json_object.find(ecut_str) != json_object.end())
    {
        if (json_object[ecut_str].is_number())
        {
            output_object[ecut_str] = json_object[ecut_str].get<double>();
        }
        else
        {
            log_fatal("The given %s option is not a number.", ecut_str.c_str());
        }
    }
    else
    {
        log_debug("No given %s . Use default %f", ecut_str.c_str(), default_ecut);
    }

    // get vcut
    if (json_object.find(vcut_str) != json_object.end())
    {
        if (json_object[vcut_str].is_number())
        {
            output_object[vcut_str] = json_object[vcut_str].get<double>();
        }
        else
        {
            log_fatal("The given %s option is not a number.", vcut_str.c_str());
        }
    }
    else
    {
        log_debug("No given %s . Use default %f", vcut_str.c_str(), default_vcut);
    }

    // get continuous randomization
    if (json_object.find(cont_str) != json_object.end())
    {
        if (json_object[cont_str].is_boolean())
        {
            output_object[cont_str] = json_object[cont_str].get<bool>();
        }
        else
        {
            log_fatal("The given %s option is not a number.", cont_str.c_str());
        }
    }
    else
    {
        log_debug("No given %s . Use default %s", cont_str.c_str(), default_contrand ? "true" : "false");
    }

    std::string output_object_str = output_object.dump();
    return output_object_str;
}


