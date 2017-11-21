/*
 * Propagator.cxx
 *
 *  Created on: 23.04.2013
 *      Author: koehne
 */

// #include <cmath>

#include <boost/property_tree/json_parser.hpp>

#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/geometry/Sphere.h"
#include "PROPOSAL/geometry/Box.h"
#include "PROPOSAL/geometry/Cylinder.h"

#include "PROPOSAL/Constants.h"

using namespace std;
using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Global defaults
// ------------------------------------------------------------------------- //

const double Propagator::global_ecut_inside_ = 500;
const double Propagator::global_ecut_infront_ = -1;
const double Propagator::global_ecut_behind_ = -1;
const double Propagator::global_vcut_inside_ = -1;
const double Propagator::global_vcut_infront_ = 0.001;
const double Propagator::global_vcut_behind_ = -1;
const double Propagator::global_cont_inside_ = false;
const double Propagator::global_cont_infront_ = true;
const double Propagator::global_cont_behind_ = false;
const bool Propagator::interpolate_ = true;  //!< Enable interpolation

// ------------------------------------------------------------------------- //
// Constructors & destructor
// ------------------------------------------------------------------------- //

Propagator::Propagator()
    : seed_(1)
    , current_sector_(NULL)
    // , particle_(MuMinusDef::Get())
    , particle_(MuMinusDef::Get())
    , detector_(new Sphere(Vector3D(), 1e18, 0))
{
    // TODO(mario): set defaults Tue 2017/09/19
    // Sector::Definition sector_def;
    // sector_def.location = Sector::ParticleLocation::InsideDetector;
    //
    // current_sector_ = new Sector(particle_);
    //
    // sectors_.push_back(current_sector_);
}

// ------------------------------------------------------------------------- //
Propagator::Propagator(const std::vector<Sector*>& sectors, const Geometry& geometry) try
    : seed_(1),
      current_sector_(NULL),
      particle_(sectors.at(0)->GetParticle()),
      detector_(geometry.clone())
{
    // --------------------------------------------------------------------- //
    // Check if all ParticleDefs are the same
    // --------------------------------------------------------------------- //

    for (std::vector<Sector*>::const_iterator iter = sectors.begin(); iter != sectors.end(); ++iter)
    {
        if ((*iter)->GetParticle().GetParticleDef() != particle_.GetParticleDef())
        {
            log_fatal("The particle definitions of the sectors must be identical for proper propagation!");
        }
        else
        {
            sectors_.push_back(new Sector(**iter));
        }
    }

    current_sector_ = sectors_.at(0);
} catch (const std::out_of_range& ex)
{
    log_fatal("No Sectors are provided for the Propagator!");
}

// ------------------------------------------------------------------------- //
Propagator::Propagator(const ParticleDef& particle_def,
                       const std::vector<SectorFactory::Definition>& sector_defs,
                       const Geometry& geometry)
    : seed_(1)
    , particle_(particle_def)
    , detector_(geometry.clone())
{
    for (std::vector<SectorFactory::Definition>::const_iterator iter = sector_defs.begin(); iter != sector_defs.end();
         ++iter)
    {
        sectors_.push_back(SectorFactory::Get().CreateSector(particle_, *iter));
    }

    try
    {
        current_sector_ = sectors_.at(0);
    }
    catch (const std::out_of_range& ex)
    {
        log_fatal("No Sectors are provided for the Propagator!");
    }
}

// ------------------------------------------------------------------------- //
Propagator::Propagator(const ParticleDef& particle_def,
                       const std::vector<SectorFactory::Definition>& sector_defs,
                       const Geometry& geometry,
                       const InterpolationDef& interpolation_def)
    : seed_(1)
    , particle_(particle_def)
    , detector_(geometry.clone())
{
    for (std::vector<SectorFactory::Definition>::const_iterator iter = sector_defs.begin(); iter != sector_defs.end();
         ++iter)
    {
        sectors_.push_back(SectorFactory::Get().CreateSector(particle_, *iter, interpolation_def));
    }

    try
    {
        current_sector_ = sectors_.at(0);
    }
    catch (const std::out_of_range& ex)
    {
        log_fatal("No Sectors are provided for the Propagator!");
    }
}

// ------------------------------------------------------------------------- //
Propagator::Propagator(const Propagator& propagator)
    : seed_(propagator.seed_)
    , sectors_(propagator.sectors_.size(), NULL)
    , current_sector_(NULL)
    , particle_(propagator.particle_)
    , detector_(propagator.detector_->clone())
{
    for(unsigned int i = 0; i < propagator.sectors_.size(); ++i)
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
    : seed_(1)
    , current_sector_(NULL)
    , particle_(particle_def)
    , detector_(NULL)
{
    double global_ecut_inside  = global_ecut_inside_;
    double global_ecut_infront = global_ecut_infront_;
    double global_ecut_behind  = global_ecut_behind_;

    double global_vcut_inside  = global_vcut_inside_;
    double global_vcut_infront = global_vcut_infront_;
    double global_vcut_behind  = global_vcut_behind_;

    bool global_cont_inside  = global_cont_inside_;
    bool global_cont_infront = global_cont_infront_;
    bool global_cont_behind  = global_cont_behind_;

    bool interpolate = interpolate_;

    // Create the json parser
    boost::property_tree::ptree pt_json;
    try
    {
        std::string expanded_config_file_path = Helper::ResolvePath(config_file);
        boost::property_tree::json_parser::read_json(expanded_config_file_path, pt_json);
    }
    catch(const boost::property_tree::ptree_error &e)
    {
        log_fatal("Unable parse \"%s\" as json file", config_file.c_str());
    }

    SetMember(seed_, "global.seed", pt_json);

    // Read in global cut and continous randomization options
    SetMember(global_ecut_inside, "global.ecut_inside", pt_json);
    SetMember(global_ecut_infront, "global.ecut_infront", pt_json);
    SetMember(global_ecut_behind, "global.ecut_behind", pt_json);
    SetMember(global_vcut_inside, "global.vcut_inside", pt_json);
    SetMember(global_vcut_infront, "global.vcut_infront", pt_json);
    SetMember(global_vcut_behind, "global.vcut_behind", pt_json);
    SetMember(global_cont_inside, "global.cont_inside", pt_json);
    SetMember(global_cont_infront, "global.cont_infront", pt_json);
    SetMember(global_cont_behind, "global.cont_behind", pt_json);

    // Read in interpolation options
    InterpolationDef interpolation_def;

    SetMember(interpolate, "global.interpolate", pt_json);
    SetMember(interpolation_def.path_to_tables, "global.path_to_tables", pt_json);
    SetMember(interpolation_def.raw, "global.raw", pt_json);

    // Read in the detector geometry
    // detector_ = GeometryFactory::Get().CreateGeometry(pt_json.get_child("detector"));

    try
    {
        boost::property_tree::ptree& detector = pt_json.get_child("detector");
        detector_ = ParseGeometryConifg(detector);
    }
    catch(const boost::property_tree::ptree_error &e)
    {
        log_fatal("You need to specify a detector geometry: %s", e.what());
    }
    // detector_ = ParseGeometryConifg(pt_json.get_child("detector"));

    // Read in global sector definition
    SectorFactory::Definition sec_def_global;

    SetMember(sec_def_global.utility_def.brems_def.multiplier, "global.brems_multiplier", pt_json);
    SetMember(sec_def_global.utility_def.photo_def.multiplier, "global.photo_multiplier", pt_json);
    SetMember(sec_def_global.utility_def.ioniz_def.multiplier, "global.ioniz_multiplier", pt_json);
    SetMember(sec_def_global.utility_def.epair_def.multiplier, "global.epair_multiplier", pt_json);
    SetMember(sec_def_global.utility_def.epair_def.lpm_effect, "global.lpm", pt_json);
    SetMember(sec_def_global.utility_def.brems_def.lpm_effect, "global.lpm", pt_json);
    SetMember(sec_def_global.do_exact_time_calculation, "global.exact_time", pt_json);
    SetMember(sec_def_global.stopping_decay, "global.stopping_decay", pt_json);

    std::string scattering = ScatteringFactory::Get().GetStringFromEnum(sec_def_global.scattering_model);
    SetMember(scattering, "global.scattering", pt_json);
    sec_def_global.scattering_model = ScatteringFactory::Get().GetEnumFromString(scattering);

    std::string brems = BremsstrahlungFactory::Get().GetStringFromEnum(sec_def_global.utility_def.brems_def.parametrization);
    SetMember(brems, "global.brems", pt_json);
    sec_def_global.utility_def.brems_def.parametrization = BremsstrahlungFactory::Get().GetEnumFromString(brems);

    std::string photo = PhotonuclearFactory::Get().GetStringFromEnum(sec_def_global.utility_def.photo_def.parametrization);
    SetMember(photo, "global.photo", pt_json);
    sec_def_global.utility_def.photo_def.parametrization = PhotonuclearFactory::Get().GetEnumFromString(photo);

    std::string shadow = PhotonuclearFactory::Get().GetStringFromShadowEnum(sec_def_global.utility_def.photo_def.shadow);
    SetMember(shadow, "global.photo_shadow", pt_json);
    sec_def_global.utility_def.photo_def.shadow = PhotonuclearFactory::Get().GetShadowEnumFromString(shadow);

    SetMember(sec_def_global.utility_def.photo_def.hardbb, "global.photo_hard_component", pt_json);

    // Read in all sector definitions
    boost::property_tree::ptree* sectors = NULL;

    try
    {
        sectors = &pt_json.get_child("sectors");
    }
    catch(const boost::property_tree::ptree_error &e)
    {
        log_fatal("You need to specify at least one sector: %s", e.what());
    }

    for (boost::property_tree::ptree::const_iterator it = sectors->begin(); it != sectors->end(); ++it)
    {
        boost::property_tree::ptree subtree = it->second;

        // Create medium
        MediumFactory::Definition medium_def;
        SetMember(medium_def.density_correction, "density_correction", subtree);

        std::string medium_name = MediumFactory::Get().GetStringFromEnum(medium_def.type);
        SetMember(medium_name, "medium", subtree);
        Medium* med = MediumFactory::Get().CreateMedium(medium_name, medium_def.density_correction);

        // Create Geometry
        Geometry* geometry = NULL;

        try
        {
            boost::property_tree::ptree& geometry_tree = subtree.get_child("geometry");
            geometry = ParseGeometryConifg(geometry_tree);
        }
        catch(const boost::property_tree::ptree_error &e)
        {
            log_fatal("You need to specify a geometry for each sector: %s", e.what());
        }

        double hirarchy = geometry->GetHirarchy();
        SetMember(hirarchy, "hirarchy", subtree);
        geometry->SetHirarchy(hirarchy);

        // Use global options in case they will not be overriden
        Sector::Definition sec_def_infront = sec_def_global;
        sec_def_infront.location = Sector::ParticleLocation::InfrontDetector;

        Sector::Definition sec_def_inside  = sec_def_global;
        sec_def_inside.location = Sector::ParticleLocation::InsideDetector;

        Sector::Definition sec_def_behind  = sec_def_global;
        sec_def_behind.location = Sector::ParticleLocation::BehindDetector;

        EnergyCutSettings cuts_infront;
        EnergyCutSettings cuts_inside;
        EnergyCutSettings cuts_behind;

        boost::optional<const boost::property_tree::ptree&> child_cuts_infront =
            it->second.get_child_optional("cuts_infront");

        if (child_cuts_infront)
        {
            cuts_infront.SetEcut(child_cuts_infront.get().get<double>("e_cut"));
            cuts_infront.SetVcut(child_cuts_infront.get().get<double>("v_cut"));
            sec_def_infront.do_continuous_randomization = child_cuts_infront.get().get<bool>("cont_rand");

        } else
        {
            cuts_infront.SetEcut(global_ecut_infront);
            cuts_infront.SetVcut(global_vcut_infront);
            sec_def_infront.do_continuous_randomization = global_cont_infront;
        }

        boost::optional<const boost::property_tree::ptree&> child_cuts_inside =
            it->second.get_child_optional("cuts_inside");

        if (child_cuts_inside)
        {
            cuts_inside.SetEcut(child_cuts_inside.get().get<double>("e_cut"));
            cuts_inside.SetVcut(child_cuts_inside.get().get<double>("v_cut"));
            sec_def_inside.do_continuous_randomization = child_cuts_inside.get().get<bool>("cont_rand");

        } else
        {
            cuts_inside.SetEcut(global_ecut_inside);
            cuts_inside.SetVcut(global_vcut_inside);
            sec_def_inside.do_continuous_randomization = global_cont_inside;
        }

        boost::optional<const boost::property_tree::ptree&> child_cuts_behind =
            it->second.get_child_optional("cuts_behind");

        if (child_cuts_behind)
        {
            cuts_behind.SetEcut(child_cuts_behind.get().get<double>("e_cut"));
            cuts_behind.SetVcut(child_cuts_behind.get().get<double>("v_cut"));
            sec_def_behind.do_continuous_randomization = child_cuts_behind.get().get<bool>("cont_rand");

        } else
        {
            cuts_behind.SetEcut(global_ecut_behind);
            cuts_behind.SetVcut(global_vcut_behind);
            sec_def_behind.do_continuous_randomization = global_cont_behind;
        }

        if (interpolate)
        {
            sectors_.push_back(new Sector(particle_, *med, cuts_infront, *geometry, sec_def_infront, interpolation_def));
            sectors_.push_back(new Sector(particle_, *med, cuts_inside, *geometry, sec_def_inside, interpolation_def));
            sectors_.push_back(new Sector(particle_, *med, cuts_behind, *geometry, sec_def_behind, interpolation_def));
        } else
        {
            sectors_.push_back(new Sector(particle_, *med, cuts_infront, *geometry, sec_def_infront));
            sectors_.push_back(new Sector(particle_, *med, cuts_inside, *geometry, sec_def_inside));
            sectors_.push_back(new Sector(particle_, *med, cuts_behind, *geometry, sec_def_behind));
        }

        delete geometry;
        delete med;
    }
}

Propagator::~Propagator()
{
    for (std::vector<Sector*>::const_iterator iter = sectors_.begin(); iter != sectors_.end(); ++iter)
    {
        delete *iter;
    }

    sectors_.clear();

    delete detector_;
}

// ------------------------------------------------------------------------- //
// Operators & swap
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
Propagator& Propagator::operator=(const Propagator &propagator)
{
    if (this != &propagator)
    {
      Propagator tmp(propagator);
      swap(tmp);
    }
    return *this;
}

// ------------------------------------------------------------------------- //
bool Propagator::operator==(const Propagator &propagator) const
{
    if (seed_ != propagator.seed_)
        return false;
    return true;
}


// ------------------------------------------------------------------------- //
bool Propagator::operator!=(const Propagator &propagator) const
{
    return !(*this == propagator);
}


// ------------------------------------------------------------------------- //
void Propagator::swap(Propagator &propagator)
{
    using std::swap;

    swap(seed_, propagator.seed_);
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

    if(Output::store_in_ASCII_file_)Output::getInstance().StorePrimaryInASCII(&particle_);

    double distance = 0;
    double result   = 0;

    // These two variables are needed to calculate the energy loss inside the detector
    // energy_at_entry_point is initialized with the current energy because this is a
    // reasonable value for particle_ which starts inside the detector

    double energy_at_entry_point = particle_.GetEnergy();
    double energy_at_exit_point  = 0;

    Vector3D particle_position = particle_.GetPosition();
    Vector3D particle_direction = particle_.GetDirection();

    bool starts_in_detector =   detector_->IsInside(particle_position, particle_direction);
    bool is_in_detector     =   false;
    bool was_in_detector    =   false;

    while(1)
    {
        particle_position = particle_.GetPosition();
        particle_direction = particle_.GetDirection();

        ChooseCurrentCollection(particle_position, particle_direction);

        if(current_sector_ == NULL)
        {
            log_info("particle_ reached the border");
            break;
        }

        // Check if have to propagate the particle_ through the whole collection
        // or only to the collection border
        distance = CalculateEffectiveDistance(particle_position, particle_direction);


        is_in_detector = detector_->IsInside(particle_position, particle_direction);
        // entry point of the detector
        if(!starts_in_detector && !was_in_detector && is_in_detector)
        {
            particle_.SetEntryPoint(particle_position);
            particle_.SetEntryEnergy( particle_.GetEnergy() );
            particle_.SetEntryTime( particle_.GetTime() );

            energy_at_entry_point = particle_.GetEnergy();

            was_in_detector = true;
        }
        // exit point of the detector
        else if(was_in_detector && !is_in_detector)
        {
            particle_.SetExitPoint(particle_position);
            particle_.SetExitEnergy( particle_.GetEnergy() );
            particle_.SetExitTime( particle_.GetTime() );

            energy_at_exit_point = particle_.GetEnergy();
            //we don't want to run in this case a second time so we set was_in_detector to false
            was_in_detector = false;

        }
        // if particle_ starts inside the detector we only ant to fill the exit point
        else if(starts_in_detector && !is_in_detector)
        {
            particle_.SetExitPoint(particle_position);
            particle_.SetExitEnergy( particle_.GetEnergy() );
            particle_.SetExitTime( particle_.GetTime() );

            energy_at_exit_point    =   particle_.GetEnergy();
            //we don't want to run in this case a second time so we set starts_in_detector to false
            starts_in_detector  =   false;

        }
        if(MaxDistance_cm <= particle_.GetPropagatedDistance() + distance)
        {
            distance = MaxDistance_cm - particle_.GetPropagatedDistance();
        }

        result  =   current_sector_->Propagate(distance);

        if(result<=0 || MaxDistance_cm <= particle_.GetPropagatedDistance()) break;
    }

    particle_.SetElost(energy_at_entry_point - energy_at_exit_point);

    #if ROOT_SUPPORT
        Output::getInstance().StorePropagatedPrimaryInTree(&particle_);
    #endif
        if(Output::store_in_ASCII_file_)Output::getInstance().StorePropagatedPrimaryInASCII(&particle_);

    //TODO(mario): undo backup Mo 2017/04/03
    // RestoreBackup_particle();

    return Output::getInstance().GetSecondarys();
}

// ------------------------------------------------------------------------- //
void Propagator::ChooseCurrentCollection(const Vector3D& particle_position, const Vector3D& particle_direction)
{
    vector<int> crossed_collections;
    crossed_collections.resize(0);

    for(unsigned int i = 0; i < sectors_.size(); i++)
    {
        // sectors_[i]->RestoreBackup_particle();

        //TODO(mario): Is that ok to delete? Tue 2017/08/08
        // if(particle_->GetType() != sectors_[i]->GetParticle()->GetType())
        // {
        //     continue;
        // }

        if(detector_->IsInfront(particle_position, particle_direction))
        {
            if(sectors_[i]->GetLocation() != 0)
            {
                continue;
            }
            else
            {
                if(sectors_[i]->GetGeometry()->IsInside(particle_position, particle_direction))
                {
                    current_sector_ = sectors_[i];
                    crossed_collections.push_back(i);
                }
                else
                {

                }
            }
        }

        else if(detector_->IsInside(particle_position, particle_direction))
        {
            if(sectors_[i]->GetLocation() != 1)
            {
                continue;
            }
            else
            {
                if(sectors_[i]->GetGeometry()->IsInside(particle_position, particle_direction))
                {
                    current_sector_ = sectors_[i];
                    crossed_collections.push_back(i);
                }
                else
                {

                }
            }

        }

        else if(detector_->IsBehind(particle_position, particle_direction))
        {
            if(sectors_[i]->GetLocation() != 2)
            {
                continue;
            }
            else
            {
                if(sectors_[i]->GetGeometry()->IsInside(particle_position, particle_direction))
                {
                    current_sector_ = sectors_[i];
                    crossed_collections.push_back(i);
                }
                //The particle reached the border of all specified collections
                else
                {

                }
            }
        }
    }

    //No process collection was found
    if(crossed_collections.size() == 0)
    {
        current_sector_ = NULL;
        log_fatal("No Cross Section was found!!!");
    }


    //Choose current collection when multiple collections are crossed!
    //
    //Choose by hirarchy of Geometry
    //If same hirarchys are available the denser one is choosen
    //If hirarchy and density are the same then the first found is taken.
    //

    for(std::vector<int>::iterator iter = crossed_collections.begin(); iter != crossed_collections.end(); ++iter)
    {
        //Current Hirachy is bigger -> Nothing to do!
        //
        if(current_sector_->GetGeometry()->GetHirarchy() >
                sectors_[*iter]->GetGeometry()->GetHirarchy() )
        {
            continue;
        }
        //Current Hirachy is equal -> Look at the density!
        //
        else if( current_sector_->GetGeometry()->GetHirarchy() ==
                 sectors_[*iter]->GetGeometry()->GetHirarchy() )
        {
            //Current Density is bigger or same -> Nothing to do!
            //

            if( current_sector_->GetMedium()->GetMassDensity() >=
                    sectors_[*iter]->GetMedium()->GetMassDensity() )
            {
                continue;
            }

            //Current Density is smaller -> Set the new collection!
            //
            else
            {
                current_sector_ =  sectors_[*iter];
            }

        }

        //Current Hirachy is smaller -> Set the new collection!
        //
        else
        {
            current_sector_ =  sectors_[*iter];
        }
    }

    //TODO(mario): Not needed anymore Thu 2017/08/24
    // if(current_sector_ != NULL)
    // {
    //     current_sector_->SetParticle(particle_);
    // }
}

// ------------------------------------------------------------------------- //
double Propagator::CalculateEffectiveDistance(const Vector3D& particle_position, const Vector3D& particle_direction)
{
    double distance_to_collection_border = 0;
    double distance_to_detector          = 0;
    double distance_to_closest_approach  = 0;

    distance_to_collection_border = current_sector_->GetGeometry()->DistanceToBorder(particle_position, particle_direction).first;
    double tmp_distance_to_border;

    for(std::vector<Sector*>::iterator iter = sectors_.begin(); iter != sectors_.end(); ++iter)
    {
        //TODO(mario): Is that ok to delete? Tue 2017/08/08
        // if (particle_.GetType() != (*iter)->Getparticle_()->GetType())
        //     continue;

        if(detector_->IsInfront(particle_position, particle_direction))
        {
            if((*iter)->GetLocation() != 0)
                continue;
            else
            {
                if((*iter)->GetGeometry()->GetHirarchy() >= current_sector_->GetGeometry()->GetHirarchy())
                {
                    tmp_distance_to_border = (*iter)->GetGeometry()->DistanceToBorder(particle_position, particle_direction).first;
                    if(tmp_distance_to_border<=0)continue;
                    distance_to_collection_border = min(tmp_distance_to_border, distance_to_collection_border);
                }
            }
        }

        else if(detector_->IsInside(particle_position, particle_direction))
        {
            if((*iter)->GetLocation() != 1)
                continue;
            else
            {
                tmp_distance_to_border = (*iter)->GetGeometry()->DistanceToBorder(particle_position, particle_direction).first;
                if(tmp_distance_to_border<=0)continue;
                distance_to_collection_border = min(tmp_distance_to_border, distance_to_collection_border);
            }

        }

        else if(detector_->IsBehind(particle_position, particle_direction))
        {
            if((*iter)->GetLocation() != 2)
                continue;
            else
            {
                if((*iter)->GetGeometry()->GetHirarchy() >= current_sector_->GetGeometry()->GetHirarchy())
                {
                    tmp_distance_to_border = (*iter)->GetGeometry()->DistanceToBorder(particle_position, particle_direction).first;
                    if(tmp_distance_to_border<=0)continue;
                    distance_to_collection_border = min(tmp_distance_to_border, distance_to_collection_border);
                }
                //The particle_ reached the border of all specified collections
                else
                {

                }
            }
        }
    }

    distance_to_detector = detector_->DistanceToBorder(particle_position, particle_direction).first;

    distance_to_closest_approach = detector_->DistanceToClosestApproach(particle_position, particle_direction);

    if(abs(distance_to_closest_approach) < GEOMETRY_PRECISION )
    {
        particle_.SetClosestApproachPoint(particle_position);
        particle_.SetClosestApproachEnergy( particle_.GetEnergy() );
        particle_.SetClosestApproachTime( particle_.GetTime() );

        distance_to_closest_approach = 0;

    }

    if(distance_to_detector > 0)
    {
        if(distance_to_closest_approach > 0)
        {
            return min(distance_to_detector, min(distance_to_collection_border, distance_to_closest_approach ));
        }
        else
        {
            return min(distance_to_detector, distance_to_collection_border);
        }
    }
    else
    {
        if(distance_to_closest_approach > 0)
        {
            return min(distance_to_closest_approach, distance_to_collection_border);
        }
        else
        {
            return distance_to_collection_border;
        }
    }
}

// ------------------------------------------------------------------------- //
// Private member functions
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
Geometry* Propagator::ParseGeometryConifg(boost::property_tree::ptree& pt)
{
    std::string origin_str = "origin";
    std::string outer_radius_str = "outer_radius";
    std::string inner_radius_str = "inner_radius";
    std::string lenght_str = "lenght";
    std::string width_str = "width";
    std::string height_str = "height";

    std::string warning_str = "Geometry %s needs to specify \"%s\" in the config file!";

    double cm_to_meter = 100.0;

    // --------------------------------------------------------------------- //
    // Get geometry from default constructor
    // --------------------------------------------------------------------- //

    Geometry* geometry = NULL;

    try
    {
        std::string name = pt.get<std::string>("shape");
        geometry = GeometryFactory::Get().CreateGeometry(name);
    }
    catch(const std::exception& e)
    {
        log_fatal(warning_str.c_str(),  "", "shape");
    }

    // --------------------------------------------------------------------- //
    // Get the position vector from the property tree
    // --------------------------------------------------------------------- //

    double x = 0;
    double y = 0;
    double z = 0;

    try
    {
        boost::property_tree::ptree child = pt.get_child(origin_str);

        if (child.size() != 3)
        {
            log_fatal("Cannot initialize Vector3D by Property Tree!");
        }

        int i = 0;
        for(boost::property_tree::ptree::const_iterator it = child.begin(); it != child.end(); ++it)
        {
            double coord = it->second.get_value<double>() * cm_to_meter;

            switch (i)
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
            ++i;
        }
    }
    catch(const std::exception& e)
    {
        log_fatal(warning_str.c_str(),  "", origin_str.c_str());
    }

    Vector3D vec(x, y, z);

    // --------------------------------------------------------------------- //
    // Check type of geometry and set members
    // --------------------------------------------------------------------- //

    if (PROPOSAL::Sphere* sphere = dynamic_cast<PROPOSAL::Sphere*>(geometry))
    {
        double radius, inner_radius;

        try
        {
            radius = pt.get<double>(outer_radius_str) * cm_to_meter;
        }
        catch(const std::exception& e)
        {
            log_fatal(warning_str.c_str(),  sphere->GetName().c_str(), outer_radius_str.c_str());
        }

        try
        {
            inner_radius = pt.get<double>(inner_radius_str) * cm_to_meter;
        }
        catch(const std::exception& e)
        {
            log_fatal(warning_str.c_str(), sphere->GetName().c_str(), inner_radius_str.c_str());
        }

        sphere->SetPosition(vec);
        sphere->SetRadius(radius);
        sphere->SetInnerRadius(inner_radius);


        return sphere;
    }
    else if (PROPOSAL::Box* box = dynamic_cast<PROPOSAL::Box*>(geometry))
    {
        double x, y, z;

        try
        {
            x = pt.get<double>(lenght_str) * cm_to_meter;
        }
        catch(const std::exception& e)
        {
            log_fatal(warning_str.c_str(), box->GetName().c_str(), lenght_str.c_str());
        }

        try
        {
            y = pt.get<double>(width_str) * cm_to_meter;
        }
        catch(const std::exception& e)
        {
            log_fatal(warning_str.c_str(), box->GetName().c_str(), width_str.c_str());
        }

        try
        {
            z = pt.get<double>(height_str) * cm_to_meter;
        }
        catch(const std::exception& e)
        {
            log_fatal(warning_str.c_str(), box->GetName().c_str(), height_str.c_str());
        }

        box->SetPosition(vec);
        box->SetX(x);
        box->SetY(y);
        box->SetZ(z);

        return box;
    }
    else if (PROPOSAL::Cylinder* cylinder = dynamic_cast<PROPOSAL::Cylinder*>(geometry))
    {
        double radius, inner_radius, z;

        try
        {
            radius = pt.get<double>(outer_radius_str) * cm_to_meter;
        }
        catch(const std::exception& e)
        {
            log_fatal(warning_str.c_str(), cylinder->GetName().c_str(), outer_radius_str.c_str());
        }

        try
        {
            inner_radius = pt.get<double>(inner_radius_str) * cm_to_meter;
        }
        catch(const std::exception& e)
        {
            log_fatal(warning_str.c_str(), cylinder->GetName().c_str(), inner_radius_str.c_str());
        }
        try
        {
            z = pt.get<double>(height_str) * cm_to_meter;
        }
        catch(const std::exception& e)
        {
            log_fatal(warning_str.c_str(), cylinder->GetName().c_str(), height_str.c_str());
        }

        cylinder->SetPosition(vec);
        cylinder->SetRadius(radius);
        cylinder->SetInnerRadius(inner_radius);
        cylinder->SetZ(z);

        return cylinder;
    }
    else
    {
        log_fatal("Dynamic casts of Geometries failed. Should not end here!");
        return NULL;
    }
}
