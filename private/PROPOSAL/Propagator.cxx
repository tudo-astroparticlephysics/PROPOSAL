/*
 * Propagator.cxx
 *
 *  Created on: 23.04.2013
 *      Author: koehne
 */

// #include <cmath>

// #include <boost/lexical_cast.hpp>
// #include <boost/algorithm/string/predicate.hpp>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "PROPOSAL/Propagator.h"

#include "PROPOSAL/sector/Sector.h"

#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"

#include "PROPOSAL/geometry/GeometryFactory.h"
#include "PROPOSAL/geometry/Sphere.h"

#include "PROPOSAL/crossection/factories/BremsstrahlungFactory.h"
#include "PROPOSAL/crossection/factories/PhotonuclearFactory.h"

#include "PROPOSAL/Output.h"
// #include "PROPOSAL/methods.h"
#include "PROPOSAL/Constants.h"
// #include "PROPOSAL/Geometry.h"

using namespace std;
using namespace PROPOSAL;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------public member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

const double Propagator::global_ecut_inside_ = 500;
const double Propagator::global_ecut_infront_ = -1;
const double Propagator::global_ecut_behind_ = -1;
const double Propagator::global_vcut_inside_ = -1;
const double Propagator::global_vcut_infront_ = 0.001;
const double Propagator::global_vcut_behind_ = -1;
const double Propagator::global_cont_inside_ = false;
const double Propagator::global_cont_infront_ = true;
const double Propagator::global_cont_behind_ = false;

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
Propagator::Propagator(const ParticleDef& particle_def, const std::string& config_file)
    : seed_(1)
    , current_sector_(NULL)
    // , particle_(particle_def)
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

    //TODO(mario): Set to global options Sun 2017/09/24
    int brems = static_cast<int>(BremsstrahlungFactory::KelnerKokoulinPetrukhin);
    int photo = static_cast<int>(PhotonuclearFactory::AbramowiczLevinLevyMaor97);
    //TODO(mario): add shadow Sun 2017/09/24
    int photo_shadow = static_cast<int>(PhotonuclearFactory::ShadowButkevichMikhailov);
    bool photo_hardbb = true;

    bool interpolate = true;

    std::string scattering = "moliere";

    // Create the json parser
    boost::property_tree::ptree pt_json;
    boost::property_tree::json_parser::read_json(config_file, pt_json);

    // Read in global options
    SetMember(seed_, "global.seed", pt_json);
    SetMember(global_ecut_inside, "global.ecut_inside", pt_json);
    SetMember(global_ecut_infront, "global.ecut_infront", pt_json);
    SetMember(global_ecut_behind, "global.ecut_behind", pt_json);
    SetMember(global_vcut_inside, "global.vcut_inside", pt_json);
    SetMember(global_vcut_infront, "global.vcut_infront", pt_json);
    SetMember(global_vcut_behind, "global.vcut_behind", pt_json);
    SetMember(global_cont_inside, "global.cont_inside", pt_json);
    SetMember(global_cont_infront, "global.cont_infront", pt_json);
    SetMember(global_cont_behind, "global.cont_behind", pt_json);

    SetMember(brems, "global.brems", pt_json);
    SetMember(photo, "global.photo", pt_json);
    SetMember(photo_shadow, "global.photo_shadow", pt_json);
    SetMember(photo_hardbb, "global.photo_hardbb", pt_json);

    SetMember(interpolate, "global.interpolate", pt_json);

    // Read in the detector geometry
    detector_ = GeometryFactory::Get().CreateGeometry(pt_json.get_child("detector"));

    // Read in global sector definition
    // SectorDef sec_def_global;
    SectorFactory::Definition sec_def_global;
    InterpolationDef interpolation_def;

    SetMember(sec_def_global.utility_def.brems_def.multiplier, "global.brems_multiplier", pt_json);
    SetMember(sec_def_global.utility_def.photo_def.multiplier, "global.photo_multiplier", pt_json);
    SetMember(sec_def_global.utility_def.ioniz_def.multiplier, "global.ioniz_multiplier", pt_json);
    SetMember(sec_def_global.utility_def.epair_def.multiplier, "global.epair_multiplier", pt_json);
    SetMember(sec_def_global.utility_def.epair_def.lpm_effect, "global.lpm", pt_json);
    SetMember(sec_def_global.utility_def.brems_def.lpm_effect, "global.lpm", pt_json);
    SetMember(sec_def_global.do_exact_time_calculation, "global.exact_time", pt_json);
    SetMember(interpolation_def.path_to_tables, "global.path_to_tables", pt_json);
    SetMember(interpolation_def.raw, "global.raw", pt_json);

    sec_def_global.scattering_model = ScatteringFactory::Get().GetEnumFromString(pt_json.get<std::string>("global.scattering"));
    sec_def_global.utility_def.brems_def.parametrization = static_cast<BremsstrahlungFactory::Enum>(brems);
    sec_def_global.utility_def.photo_def.parametrization = static_cast<PhotonuclearFactory::Enum>(photo);
    sec_def_global.utility_def.photo_def.shadow = static_cast<PhotonuclearFactory::Shadow>(photo_shadow);
    sec_def_global.utility_def.photo_def.hardbb = photo_hardbb;

    // Read in all sector definitions
    boost::property_tree::ptree sectors = pt_json.get_child("sectors");

    for (boost::property_tree::ptree::const_iterator it = sectors.begin(); it != sectors.end(); ++it)
    {
        Medium* med        = MediumFactory::Get().CreateMedium(it->second.get<std::string>("medium"));

        Geometry* geometry = GeometryFactory::Get().CreateGeometry(it->second.get_child("geometry"));
        geometry->SetHirarchy(it->second.get<unsigned int>("hirarchy"));

        // Use global options in case they will not be overriden
        SectorFactory::Definition sec_def_infront = sec_def_global;
        sec_def_infront.location = Sector::ParticleLocation::InfrontDetector;

        SectorFactory::Definition sec_def_inside  = sec_def_global;
        sec_def_inside.location = Sector::ParticleLocation::InsideDetector;

        SectorFactory::Definition sec_def_behind  = sec_def_global;
        sec_def_behind.location = Sector::ParticleLocation::BehindDetector;

        double density_correction = it->second.get<double>("density_correction");

        sec_def_infront.medium_def.density_correction = density_correction;
        sec_def_inside.medium_def.density_correction  = density_correction;
        sec_def_behind.medium_def.density_correction  = density_correction;

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

void Propagator::ChooseCurrentCollection(Vector3D& particle_position, Vector3D& particle_direction)
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


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

double Propagator::CalculateEffectiveDistance(Vector3D& particle_position, Vector3D& particle_direction)
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


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// double Propagator::Propagate( double distance )
// {
//     bool    flag;
//     double  displacement;
//
//     double propagated_distance  =   0;
//
//     double  initial_energy  =   particle_->GetEnergy();
//     double  final_energy    =   particle_->GetEnergy();
//
//     pair<double, ParticleType::Enum> decay;
//     std::vector<PROPOSALParticle*> decay_products;
//
//     pair<double, ParticleType::Enum> energy_loss;
//
//
//     int secondary_id    =   0;
//
//     //first: final energy befor first interaction second: energy at which the
//     // particle decay
//     //first and second are compared to decide if interaction happens or decay
//     pair<double,double> energy_till_stochastic_;
//
//
//     if(distance < 0)
//     {
//         distance   =   0;
//     }
//
//     if(initial_energy <= particle_->GetLow() || distance==0)
//     {
//         flag    =   false;
//     }
//     else
//     {
//         flag    =   true;
//     }
//
//     while(flag)
//     {
//         energy_till_stochastic_ = current_sector_->CalculateEnergyTillStochastic(particle_, initial_energy );
//         if(energy_till_stochastic_.first > energy_till_stochastic_.second)
//         {
//             particle_interaction_   =   true;
//             final_energy            =   energy_till_stochastic_.first;
//         }
//         else
//         {
//             particle_interaction_   =   false;
//             final_energy            =   energy_till_stochastic_.second;
//
//         }
//
//         //Calculate the displacement according to initial energy initial_energy and final_energy
//         displacement  =   current_sector_->CalculateDisplacement(
//                     initial_energy,
//                     final_energy,
//                     current_sector_->GetDensityCorrection()*(distance - propagated_distance)) / current_sector_->GetDensityCorrection();
//
//         // The first interaction or decay happens behind the distance we want to propagate
//         // So we calculate the final energy using only continuous losses
//         if( displacement > distance - propagated_distance )
//         {
//             displacement  =   distance - propagated_distance;
//
//             final_energy  =   current_sector_->CalculateFinalEnergy(particle_, initial_energy, current_sector_->GetDensityCorrection()*displacement);
//
//         }
//         //Advance the Particle according to the displacement
//         //Initial energy and final energy are needed if Molier Scattering is enabled
//         AdvanceParticle(displacement, initial_energy, final_energy);
//
//         propagated_distance +=  displacement;
//
//         if(abs(distance - propagated_distance) < abs(distance)*COMPUTER_PRECISION)
//         {
//             propagated_distance = distance;  // computer precision control
//         }
//         //Randomize the continuous energy loss if this option is enabled
//         if( current_sector_->GetDoRandomization() )
//         {
//             if(final_energy != particle_->GetLow())
//             {
//                 final_energy  = current_sector_->Randomize( initial_energy, final_energy );
//             }
//
//         }
//         // Lower limit of particle energy is reached or
//         // or complete particle is propagated the whole distance
//         if( final_energy == particle_->GetLow() || propagated_distance == distance)
//         {
//             break;
//         }
//
//         //Set the particle energy to the current energy before making
//         //stochatic losses or decay
//         particle_->SetEnergy( final_energy );
//
//         if(particle_interaction_)
//         {
//             energy_loss     =   current_sector_->MakeStochasticLoss(particle_);
//             if (energy_loss.second == ParticleType::unknown)
//             {
//                 // in this case, no cross section is chosen, so there is no interaction
//                 // due to the parameterization of the cross section cutoffs
//                 log_debug("no interaction due to the parameterization of the cross section cutoffs. final energy: %f\n", final_energy);
//                 initial_energy = final_energy;
//                 continue;
//             }
//             final_energy    -=  energy_loss.first;
//             // log_debug("Energyloss: %f\t%s", energy_loss.first, PROPOSALParticle::GetName(energy_loss.second).c_str());
//             secondary_id    =   particle_->GetParticleId() + 1;
//             Output::getInstance().FillSecondaryVector(particle_, secondary_id, energy_loss, 0);
//         }
//         else
//         {
//             DecayChannel* mode = &particle_->GetDecayTable().SelectChannel();
//             decay_products = mode->Decay(particle_);
//             Output::getInstance().FillSecondaryVector(decay_products);
//
//             //TODO(mario): Delete decay products Tue 2017/08/22
//
//             // decay           =   current_sector_->MakeDecay();
//             // final_energy    =   0;
//             // log_debug("Decay of particle: %s", particle_->GetName().c_str());
//             // secondary_id    = particle_->GetParticleId()  +   1;
//             // Output::getInstance().FillSecondaryVector(particle_, secondary_id, decay ,0);
//
//         }
//
//         //break if the lower limit of particle energy is reached
//         if(final_energy <= particle_->GetLow())
//         {
//             break;
//         }
//
//         //Next round: update the inital energy
//         initial_energy  =   final_energy;
//
//     }
//
//     // if(stopping_decay_)
//     // {
//     //     if(propagated_distance!=distance && final_energy!=0 && particle_->GetLifetime()>=0)
//     //     {
//     //         particle_->SetEnergy(particle_->GetMass());
//     //
//     //         double t    =   particle_->GetT() -particle_->GetLifetime()*log(RandomDouble());
//     //         double product_energy   =   0;
//     //
//     //         pair<double, ParticleType::Enum> decay_to_store;
//     //         secondary_id    =   particle_->GetParticleId() + 1;
//     //
//     //         particle_->SetT( t );
//     //
//     //         if(particle_->GetType()==2)
//     //         {
//     //             // --------------------------------------------------------------------- //
//     //             // Calculate random numbers before passing to a fuction, because
//     //             // the order of argument evaluation is unspecified in c++ standards and
//     //             // therfor depend on the compiler.
//     //             // --------------------------------------------------------------------- //
//     //
//     //             double rnd1 = RandomDouble();
//     //             double rnd2 = RandomDouble();
//     //
//     //             product_energy  =   current_sector_->GetDecay()->CalculateProductEnergy(rnd1, 0.5, rnd2);
//     //         }
//     //         else
//     //         {
//     //             product_energy  =   current_sector_->GetDecay()->CalculateProductEnergy(RandomDouble(), 0.5, 0.5);
//     //         }
//     //
//     //         decay_to_store.first    =   product_energy;
//     //         decay_to_store.second   =   current_sector_->GetDecay()->GetOut();
//     //
//     //         final_energy  =   0;
//     //
//     //         Output::getInstance().FillSecondaryVector(particle_,secondary_id, decay_to_store, final_energy);
//     //     }
//     // }
//
//     particle_->SetEnergy(final_energy);
//
//     //Particle reached the border, final energy is returned
//     if(propagated_distance==distance)
//     {
//         return final_energy;
//     }
//     //The particle stopped/decayed, the propageted distance is return with a minus sign
//     else
//     {
//         return -propagated_distance;
//     }
//     //Should never be here
//     return 0;
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// void Propagator::AdvanceParticle(double dr, double ei, double ef)
// {
//
//     double dist = particle_->GetPropagatedDistance();
//     double time = particle_->GetT();
//     Vector3D position = particle_->GetPosition();
//
//     dist   +=  dr;
//
//     if(do_exact_time_calculation_)
//     {
//         time   +=  current_sector_->CalculateParticleTime(ei, ef)/current_sector_->GetDensityCorrection();
//     }
//     else
//     {
//         time   +=  dr/SPEED;
//     }
//
//
//     if(scattering_model_!=-1)
//     {
//         switch(scattering_model_)
//         {
//             case 0:
//                 current_sector_->GetScattering()->Scatter(dr,ei,ef);
//                 break;
//
//             case 1:
//                 scatteringFirstOrder_->Scatter(dr, particle_, current_sector_->GetMedium());
//                 break;
//
//             case 2:
//                 scatteringFirstOrderMoliere_->Scatter(dr, particle_, current_sector_->GetMedium());
//                 break;
//             default:
//                 log_error("Never should be here! scattering_model = %i !",scattering_model_);
//         }
//
//     }
//     else
//     {
//         position = position + dr*particle_->GetDirection();
//         particle_->SetPosition(position);
//     }
//     particle_->SetPropagatedDistance(dist);
//     particle_->SetT(time);
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// void Propagator::ReadConfigFile(string config_file, bool DoApplyOptions)
// {
//     bool found_detector         =   false;
//
//     //global
//
//     if(!FileExist(config_file))
//     {
//         log_fatal("Error: config file %s does not exist!",config_file.c_str());
//         exit(1);
//     }
//
//     ifstream file;
//     file.open(config_file.c_str());
//
//    // int mediacur;
//     char buf[256];
//     string str_buf;
//     deque<string> *token;
//     string taux;
//
//     while(file.good())
//     {
//         file.getline(buf,256);
//
//         str_buf =   string(buf);
//         token   =   SplitString(str_buf," \t");
//
//         // Ignore empty lines
//         if(token->empty())
//         {
//             continue;
//         }
//
//         taux =   NextToken(token);
//
//         // Ignore lines starting with # (comments)
//         if(taux.at(0)=='#')
//         {
//             continue;
//         }
//
//         // Reading the global options
//         else if(ToLowerCase(taux).compare("global")==0)
//         {
//             log_info("Reading the global options");
//             continue;
//         }
//         // seed
//         else if(ToLowerCase(taux).compare("seed")==0)
//         {
//             taux    =   NextToken(token);
//
//             try {
//                 seed_ = boost::lexical_cast<int>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("The seed is set to %s but must be an integer! Set to 1", taux.c_str());
//                 seed_ = 1;
//             }
//         }
//         // bremsstrahlungs parametrization
//         else if(ToLowerCase(taux).compare("brems")==0)
//         {
//             taux    =   NextToken(token);
//
//             try {
//                 brems_ = static_cast<ParametrizationType::Enum>(boost::lexical_cast<int>(taux));
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("The bremsstrahlungs parametrization indentifier is set to %s \
//                     but must be a ParametrizationType::Enum! Set to BremsKelnerKokoulinPetrukhin", taux.c_str());
//                 brems_ = ParametrizationType::BremsKelnerKokoulinPetrukhin;
//             }
//         }
//         // photonuclear parametrization
//         else if(ToLowerCase(taux).compare("photo")==0)
//         {
//             taux    =   NextToken(token);
//
//             try {
//                 photo_ = static_cast<ParametrizationType::Enum>(boost::lexical_cast<int>(taux));
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("The photonuclear parametrization indentifier is set to %s \
//                     but must be a ParametrizationType::Enum! Set to PhotoAbramowiczLevinLevyMaor97ShadowButkevich", taux.c_str());
//                 photo_ = ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich;
//             }
//         }
//         // bremsstahlungs mulitpiler
//         else if(ToLowerCase(taux).compare("brems_multiplier")==0)
//         {
//             taux    =   NextToken(token);
//
//             try {
//                 brems_multiplier_  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("The bremsstrahlungs multiplier is set to %s but must be a double! Set to 1.", taux.c_str());
//                 brems_multiplier_ = 1.;
//             }
//         }
//         // photonuclear multiplier
//         else if(ToLowerCase(taux).compare("photo_multiplier")==0)
//         {
//             taux    =   NextToken(token);
//
//             try {
//                 photo_multiplier_  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("The photonuclear multiplier is set to %s but must be a double! Set to 1.", taux.c_str());
//                 photo_multiplier_ = 1.;
//             }
//         }
//         // epairproduction multiplier
//         else if(ToLowerCase(taux).compare("epair_multiplier")==0)
//         {
//             taux    =   NextToken(token);
//
//             try {
//                 epair_multiplier_  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("The epairproduction multiplier is set to %s but must be a double! Set to 1.", taux.c_str());
//                 epair_multiplier_ = 1.;
//             }
//         }
//         // ionization multiplier
//         else if(ToLowerCase(taux).compare("ioniz_multiplier")==0)
//         {
//             taux    =   NextToken(token);
//
//             try {
//                 ioniz_multiplier_  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("The ionization multiplier is set to %s but must be a double! Set to 1.", taux.c_str());
//                 ioniz_multiplier_ = 1.;
//             }
//         }
//         // global ecut inside the detector
//         else if(ToLowerCase(taux).compare("ecut_inside")==0)
//         {
//             taux    =   NextToken(token);
//
//             try {
//                 global_ecut_inside_  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("The ecut for inside the detector is set to %s but must be a double! Set to 500.", taux.c_str());
//                 global_ecut_inside_ = 500.;
//             }
//         }
//         // global ecut behind the detector
//         else if(ToLowerCase(taux).compare("ecut_behind")==0)
//         {
//             taux    =   NextToken(token);
//
//             try {
//                 global_ecut_behind_  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("The ecut for behind the detector is set to %s but must be a double! Set to -1.", taux.c_str());
//                 global_ecut_behind_ = -1.;
//             }
//         }
//         // global ecut infront of the detector
//         else if(ToLowerCase(taux).compare("ecut_infront")==0)
//         {
//             taux    =   NextToken(token);
//
//             try {
//                 global_ecut_infront_  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("The ecut for infront of the detector is set to %s but must be a double! Set to -1.", taux.c_str());
//                 global_ecut_infront_ = -1.;
//             }
//         }
//         // global vcut inside the detector
//         else if(ToLowerCase(taux).compare("vcut_inside")==0)
//         {
//             taux    =   NextToken(token);
//
//             try {
//                 global_vcut_inside_  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("The vcut for inside the detector is set to %s but must be a double! Set to -1.", taux.c_str());
//                 global_vcut_inside_ = -1.;
//             }
//         }
//         // global vcut behind the detector
//         else if(ToLowerCase(taux).compare("vcut_behind")==0)
//         {
//             taux    =   NextToken(token);
//
//             try {
//                 global_vcut_behind_  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("The vcut for behind the detector is set to %s but must be a double! Set to -1.", taux.c_str());
//                 global_vcut_behind_ = -1.;
//             }
//         }
//         // global vcut infront of the detector
//         else if(ToLowerCase(taux).compare("vcut_infront")==0)
//         {
//             taux    =   NextToken(token);
//
//             try {
//                 global_vcut_infront_  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("The vcut for infront of the detector is set to %s but must be a double! Set to 0.001", taux.c_str());
//                 global_vcut_infront_ = 0.001;
//             }
//         }
//         // global continuous randominzation inside the detector
//         else if(ToLowerCase(taux).compare("cont_inside")==0)
//         {
//             taux    =   NextToken(token);
//
//             try {
//                 global_cont_inside_  = boost::lexical_cast<bool>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("cont for inside the detector is set to %s but must be a bool! Set to false", taux.c_str());
//                 global_cont_inside_ = false;
//             }
//         }
//         // global continuous randominzation behind the detector
//         else if(ToLowerCase(taux).compare("cont_behind")==0)
//         {
//             taux    =   NextToken(token);
//
//             try {
//                 global_cont_behind_  = boost::lexical_cast<bool>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("cont for behind the detector is set to %s but must be a double! Set to false.", taux.c_str());
//                 global_vcut_behind_ = false;
//             }
//         }
//         // global continuous randominzation infront of the detector
//         else if(ToLowerCase(taux).compare("cont_infront")==0)
//         {
//             taux    =   NextToken(token);
//
//             try {
//                 global_cont_infront_  = boost::lexical_cast<bool>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("cont for infront of the detector is set to %s but must be a double! Set to true", taux.c_str());
//                 global_vcut_infront_ = true;
//             }
//         }
//         // lpm effect
//         else if(ToLowerCase(taux).compare("lpm")==0)
//         {
//             lpm_ =   true;
//         }
//         // moliere scattering
//         else if(ToLowerCase(taux).compare("moliere")==0 || ToLowerCase(taux).compare("scattering")==0)
//         {
//             if(ToLowerCase(taux).compare("moliere")==0)
//             {
//                 moliere_ =   true;
//                 scattering_model_ = 0;
//                 continue;
//             }
//             taux = NextToken(token);
//             if(ToLowerCase(taux).compare("moliere")==0)
//             {
//                 moliere_ =   true;
//                 scattering_model_ = 0;
//                 continue;
//             }
//             if(ToLowerCase(taux).compare("firstorder")==0)
//             {
//                 scatteringFirstOrder_ = new ScatteringFirstOrder();
//                 scattering_model_ = 1;
//                 continue;
//             }
//             if(ToLowerCase(taux).compare("firstordermoliere")==0)
//             {
//                 scatteringFirstOrderMoliere_ = new ScatteringMoliere();
//                 scattering_model_ = 2;
//                 continue;
//             }
//             log_error("Scattering model not known. Defaulting to moliere!");
//             moliere_ =   true;
//             scattering_model_ = 0;
//             continue;
//         }
//         // exact location time
//         else if(ToLowerCase(taux).compare("exact_time")==0)
//         {
//             do_exact_time_calculation_ =   true;
//         }
//         // do not interpolate: intergrate everything
//         else if(ToLowerCase(taux).compare("integrate")==0)
//         {
//             integrate_ =   true;
//         }
//         // save interpolation tables binary or not
//         else if(ToLowerCase(taux).compare("raw")==0)
//         {
//             raw_ =   true;
//         }
//         // path to interpolation tables
//         else if(ToLowerCase(taux).compare("path_to_tables")==0)
//         {
//             taux    =   NextToken(token);
//             path_to_tables_ =   taux;
//         }
//
//         //Builing the detector geometry
//         else if(ToLowerCase(taux).compare("detector")==0)
//         {
//             if(found_detector)
//             {
//                 log_warn("Detector already specified. There can be only one. This one will be ignored");
//                 continue;
//             }
//             found_detector  =   true;
//
//             //find the first not empty line not starting with #
//             while(file.good())
//             {
//                 file.getline(buf,256);
//                 str_buf =   string(buf);
//                 token   =   SplitString(str_buf," \t");
//                 if(!token->empty())
//                 {
//                     taux =   NextToken(token);
//                     if(taux.at(0)!='#')
//                     {
//                         break;
//                     }
//                 }
//             }
//
//             // detector_   =   new Geometry();
//             detector_ = InitGeometry(token,taux);
//
//         }
//
//         // a sector consists of geometry, a medium and cut settings
//         // here the ProccessCollections are initialized
//         else if(ToLowerCase(taux).compare("sector")==0)
//         {
//             InitProcessCollections(file);
//         }
//         else
//         {
//             log_warn("Unrecognized option: %s",taux.c_str());
//             continue;
//         }
//     }
//
//     if(DoApplyOptions)
//     {
//         ApplyOptions();
//     }
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

// Propagator::Propagator(
//                        ParticleDef particle_def,
//                        std::string path_to_tables,
//                        bool exact_time,
//                        bool lpm,
//                        bool integrate,
//                        int scattering_model)
//     :order_of_interpolation_    ( 5 )
//     ,debug_                     ( false )
//     ,seed_                      ( 1 )
//     ,brems_                     ( ParametrizationType::BremsKelnerKokoulinPetrukhin )
//     ,photo_                     ( ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich )
//     ,lpm_                       ( lpm )
//     ,stopping_decay_            ( true )
//     ,do_exact_time_calculation_ ( exact_time )
//     ,integrate_                 ( integrate )
//     ,brems_multiplier_          ( 1 )
//     ,photo_multiplier_          ( 1 )
//     ,ioniz_multiplier_          ( 1 )
//     ,epair_multiplier_          ( 1 )
//     ,global_ecut_inside_        ( 500 )
//     ,global_ecut_infront_       ( -1 )
//     ,global_ecut_behind_        ( -1 )
//     ,global_vcut_inside_        ( -1 )
//     ,global_vcut_infront_       ( 0.001 )
//     ,global_vcut_behind_        ( -1 )
//     ,global_cont_inside_        ( false )
//     ,global_cont_infront_       ( true )
//     ,global_cont_behind_        ( false )
//     ,path_to_tables_            ( path_to_tables )
//     ,raw_                       ( false )
//     ,scattering_model_          ( scattering_model )
//     ,current_sector_        (NULL)
// {
//     particle_              = new PROPOSALParticle(particle_def);
//     backup_particle_       = new PROPOSALParticle(*particle_);
//     detector_              = new Sphere(Vector3D(), 1e18, 0);
//     // detector_              = new Geometry();
//     // detector_->InitSphere(0,0,0,1e18,0);
//
//     if(scattering_model_ != -1)
//     {
//         switch(scattering_model_)
//         {
//             case 0:
//                 moliere_ = true;
//                 scattering_model_ = 0;
//                 break;
//             case 1:
//                 moliere_ = false;
//                 scatteringFirstOrder_ =   new ScatteringFirstOrder();
//                 break;
//             case 2:
//                 moliere_ = false;
//                 scatteringFirstOrderMoliere_ =   new ScatteringMoliere();
//                 break;
//             default:
//                 log_error("scattering model not known! defaulting to moliere!");
//                 moliere_ = true;
//                 scattering_model_ = 0;
//         }
//     }
// }
//
// //----------------------------------------------------------------------------//
// //----------------------------------------------------------------------------//
//
//
// Propagator::Propagator(Medium* medium,
//                        EnergyCutSettings* cuts,
//                        ParticleDef particle_def,
//                        string path_to_tables,
//                        bool moliere,
//                        bool continuous_rand,
//                        bool exact_time,
//                        bool lpm,
//                        ParametrizationType::Enum brems,
//                        ParametrizationType::Enum photo,
//                        double brems_multiplier,
//                        double photo_multiplier,
//                        double ioniz_multiplier,
//                        double epair_multiplier,
//                        bool integrate,
//                        int scattering_model
//                        )
//     :order_of_interpolation_    ( 5 )
//     ,debug_                     ( false )
//     ,particle_interaction_      ( false )
//     ,seed_                      ( 1 )
//     ,brems_                     ( brems )
//     ,photo_                     ( photo )
//     ,lpm_                       ( lpm )
//     ,moliere_                   ( moliere )
//     ,stopping_decay_            ( true )
//     ,do_exact_time_calculation_ ( exact_time )
//     ,integrate_                 ( integrate )
//     ,brems_multiplier_          ( brems_multiplier )
//     ,photo_multiplier_          ( photo_multiplier )
//     ,ioniz_multiplier_          ( ioniz_multiplier )
//     ,epair_multiplier_          ( epair_multiplier )
//     ,global_ecut_inside_        ( 500 )
//     ,global_ecut_infront_       ( -1 )
//     ,global_ecut_behind_        ( -1 )
//     ,global_vcut_inside_        ( -1 )
//     ,global_vcut_infront_       ( 0.001 )
//     ,global_vcut_behind_        ( -1 )
//     ,global_cont_inside_        ( false )
//     ,global_cont_infront_       ( true )
//     ,global_cont_behind_        ( false )
//     ,path_to_tables_            ( path_to_tables )
//     ,raw_                       ( true )
//     ,scattering_model_          (scattering_model)
//     ,current_sector_        (NULL)
// {
//     particle_              = new PROPOSALParticle(particle_def);
//     backup_particle_       = new PROPOSALParticle(*particle_);
//     current_sector_    = new ProcessCollection(particle_, medium, cuts);
//     detector_              = new Sphere(Vector3D(), 1e18, 0);
//     // detector_              = new Geometry();
//     // detector_->InitSphere(0,0,0,1e18,0);
//     current_sector_->SetGeometry(detector_);
//
//     current_sector_->SetLocation(1); // Inside the detector
//     sectors_.push_back(current_sector_);
//
//     if(continuous_rand)
//     {
//         current_sector_->EnableContinuousRandomization();
//     }
//
//     if(scattering_model_ != -1)
//     {
//         switch(scattering_model_)
//         {
//             case 0:
//                 moliere_ = true;
//                 scattering_model_ = 0;
//                 break;
//             case 1:
//                 moliere_ = false;
//                 scatteringFirstOrder_ =   new ScatteringFirstOrder();
//                 break;
//             case 2:
//                 moliere_ = false;
//                 scatteringFirstOrderMoliere_ =   new ScatteringMoliere();
//                 break;
//             default:
//                 log_error("scattering model not known! defaulting to moliere!");
//                 moliere_ = true;
//                 scattering_model_ = 0;
//         }
//     }
//
//     ApplyOptions();
// }
// //----------------------------------------------------------------------------//
// //----------------------------------------------------------------------------//
//
//
// Propagator::Propagator(string config_file, bool DoApplyOptions)
//     :order_of_interpolation_    ( 5 )
//     ,debug_                     ( false )
//     ,particle_interaction_      ( false )
//     ,seed_                      ( 1 )
//     ,brems_                     ( ParametrizationType::BremsKelnerKokoulinPetrukhin )
//     ,photo_                     ( ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich )
//     ,lpm_                       ( false )
//     ,moliere_                   ( false )
//     ,stopping_decay_            ( true )
//     ,do_exact_time_calculation_ ( false )
//     ,integrate_                 ( false )
//     ,brems_multiplier_          ( 1 )
//     ,photo_multiplier_          ( 1 )
//     ,ioniz_multiplier_          ( 1 )
//     ,epair_multiplier_          ( 1 )
//     ,global_ecut_inside_        ( 500 )
//     ,global_ecut_infront_       ( -1 )
//     ,global_ecut_behind_        ( -1 )
//     ,global_vcut_inside_        ( -1 )
//     ,global_vcut_infront_       ( 0.001 )
//     ,global_vcut_behind_        ( -1 )
//     ,global_cont_inside_        ( false )
//     ,global_cont_infront_       ( true )
//     ,global_cont_behind_        ( false )
//     ,path_to_tables_            ( "" )
//     ,raw_                       ( false )
//     ,particle_                  (NULL)
//     ,backup_particle_           (NULL)
//     ,scattering_model_          (-1)
//     ,current_sector_        (NULL)
// {
//     ReadConfigFile(config_file, DoApplyOptions);
// }
//
// //----------------------------------------------------------------------------//
// //----------------------------------------------------------------------------//
//
// Propagator::Propagator(std::string config_file, PROPOSALParticle* particle, bool DoApplyOptions)
//     :order_of_interpolation_    ( 5 )
//     ,debug_                     ( false )
//     ,particle_interaction_      ( false )
//     ,seed_                      ( 1 )
//     ,brems_                     ( ParametrizationType::BremsKelnerKokoulinPetrukhin )
//     ,photo_                     ( ParametrizationType::PhotoAbramowiczLevinLevyMaor97ShadowButkevich )
//     ,lpm_                       ( false )
//     ,moliere_                   ( false )
//     ,stopping_decay_            ( true )
//     ,do_exact_time_calculation_ ( false )
//     ,integrate_                 ( false )
//     ,brems_multiplier_          ( 1 )
//     ,photo_multiplier_          ( 1 )
//     ,ioniz_multiplier_          ( 1 )
//     ,epair_multiplier_          ( 1 )
//     ,global_ecut_inside_        ( 500 )
//     ,global_ecut_infront_       ( -1 )
//     ,global_ecut_behind_        ( -1 )
//     ,global_vcut_inside_        ( -1 )
//     ,global_vcut_infront_       ( 0.001 )
//     ,global_vcut_behind_        ( -1 )
//     ,global_cont_inside_        ( false )
//     ,global_cont_infront_       ( true )
//     ,global_cont_behind_        ( false )
//     ,path_to_tables_            ( "" )
//     ,raw_                       ( false )
//     ,scattering_model_          (-1)
//     ,current_sector_        (NULL)
// {
//     particle_        = new PROPOSALParticle(*particle);
//     backup_particle_ = new PROPOSALParticle(*particle_);
//     ReadConfigFile(config_file, DoApplyOptions);
// }

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Propagator::Propagator(const Propagator &propagator)
    :seed_                      ( propagator.seed_ )
    ,particle_                      ( propagator.particle_ )
    // ,brems_                     ( propagator.brems_ )
    // ,photo_                     ( propagator.photo_ )
    // ,lpm_                       ( propagator.lpm_ )
    // ,moliere_                   ( propagator.moliere_ )
    // ,stopping_decay_            ( propagator.stopping_decay_ )
    // ,do_exact_time_calculation_ ( propagator.do_exact_time_calculation_ )
    // ,integrate_                 ( propagator.integrate_ )
    // ,brems_multiplier_          ( propagator.brems_multiplier_ )
    // ,photo_multiplier_          ( propagator.photo_multiplier_ )
    // ,ioniz_multiplier_          ( propagator.ioniz_multiplier_ )
    // ,epair_multiplier_          ( propagator.epair_multiplier_ )
    // ,global_ecut_inside_        ( propagator.global_ecut_inside_ )
    // ,global_ecut_infront_       ( propagator.global_ecut_infront_ )
    // ,global_ecut_behind_        ( propagator.global_ecut_behind_ )
    // ,global_vcut_inside_        ( propagator.global_vcut_inside_ )
    // ,global_vcut_infront_       ( propagator.global_vcut_infront_ )
    // ,global_vcut_behind_        ( propagator.global_vcut_behind_ )
    // ,global_cont_inside_        ( propagator.global_cont_inside_ )
    // ,global_cont_infront_       ( propagator.global_cont_infront_ )
    // ,global_cont_behind_        ( propagator.global_cont_behind_ )
    // ,path_to_tables_            ( propagator.path_to_tables_ )
    // ,raw_                       ( propagator.raw_ )
    // ,particle_                  ( propagator.particle_ )
    // ,backup_particle_           ( propagator.backup_particle_ )
    //FirstOrderScattering
    // ,scatteringFirstOrder_          ( propagator.scatteringFirstOrder_ )
    // ,scatteringFirstOrderMoliere_   ( propagator.scatteringFirstOrderMoliere_ )
    // ,scattering_model_          (propagator.scattering_model_)
    // ,current_sector_        ( new Collection(*propagator.current_sector_) )

{

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


Propagator& Propagator::operator=(const Propagator &propagator)
{
    if (this != &propagator)
    {
      Propagator tmp(propagator);
      swap(tmp);
    }
    return *this;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Propagator::operator==(const Propagator &propagator) const
{
    // if( particle_                 != propagator.particle_ )               return false;
    //FirstOrderScattering
    // if( scatteringFirstOrder_           != propagator.scatteringFirstOrder_ )           return false;
    // if( scatteringFirstOrderMoliere_    != propagator.scatteringFirstOrderMoliere_ )    return false;
    // if( scattering_model_               != propagator.scattering_model_)                return false;

    if( seed_                     != propagator.seed_ )                   return false;
    // if( brems_                    != propagator.brems_ )                  return false;
    // if( photo_                    != propagator.photo_ )                  return false;
    // if( lpm_                      != propagator.lpm_ )                    return false;
    // if( moliere_                  != propagator.moliere_ )                return false;
    // if( stopping_decay_           != propagator.stopping_decay_ )         return false;
    // if( do_exact_time_calculation_!= propagator.do_exact_time_calculation_ )return false;
    // if( integrate_                != propagator.integrate_ )              return false;
    // if( brems_multiplier_         != propagator.brems_multiplier_ )       return false;
    // if( photo_multiplier_         != propagator.photo_multiplier_ )       return false;
    // if( ioniz_multiplier_         != propagator.ioniz_multiplier_ )       return false;
    // if( epair_multiplier_         != propagator.epair_multiplier_ )       return false;
    // if( global_ecut_inside_       != propagator.global_ecut_inside_ )     return false;
    // if( global_ecut_infront_      != propagator.global_ecut_infront_ )    return false;
    // if( global_ecut_behind_       != propagator.global_ecut_behind_ )     return false;
    // if( global_vcut_inside_       != propagator.global_vcut_inside_ )     return false;
    // if( global_vcut_infront_      != propagator.global_vcut_infront_ )    return false;
    // if( global_vcut_behind_       != propagator.global_vcut_behind_ )     return false;
    // if( global_cont_inside_       != propagator.global_cont_inside_ )     return false;
    // if( global_cont_infront_      != propagator.global_cont_infront_ )    return false;
    // if( global_cont_behind_       != propagator.global_cont_behind_ )     return false;
    // if( *current_sector_      != *propagator.current_sector_ )    return false;
    // if( raw_                      != propagator.raw_ )                    return false;
    //
    // if( path_to_tables_.compare( propagator.path_to_tables_ )!=0 )        return false;

    //else
    return true;
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


bool Propagator::operator!=(const Propagator &propagator) const
{
    return !(*this == propagator);
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


void Propagator::swap(Propagator &propagator)
{
    using std::swap;

    swap( seed_                     ,   propagator.seed_ );
    // swap( brems_                    ,   propagator.brems_ );
    // swap( photo_                    ,   propagator.photo_ );
    // swap( lpm_                      ,   propagator.lpm_ );
    // swap( moliere_                  ,   propagator.moliere_ );
    // swap( stopping_decay_           ,   propagator.stopping_decay_ );
    // swap( do_exact_time_calculation_,   propagator.do_exact_time_calculation_ );
    // swap( integrate_                ,   propagator.integrate_ );
    // swap( brems_multiplier_         ,   propagator.brems_multiplier_ );
    // swap( photo_multiplier_         ,   propagator.photo_multiplier_ );
    // swap( ioniz_multiplier_         ,   propagator.ioniz_multiplier_ );
    // swap( epair_multiplier_         ,   propagator.epair_multiplier_ );
    // swap( global_ecut_inside_       ,   propagator.global_ecut_inside_ );
    // swap( global_ecut_infront_      ,   propagator.global_ecut_infront_ );
    // swap( global_ecut_behind_       ,   propagator.global_ecut_behind_ );
    // swap( global_vcut_inside_       ,   propagator.global_vcut_inside_ );
    // swap( global_vcut_infront_      ,   propagator.global_vcut_infront_ );
    // swap( global_vcut_behind_       ,   propagator.global_vcut_behind_ );
    // swap( global_cont_inside_       ,   propagator.global_cont_inside_ );
    // swap( global_cont_infront_      ,   propagator.global_cont_infront_ );
    // swap( global_cont_behind_       ,   propagator.global_cont_behind_ );
    // swap( raw_                      ,   propagator.raw_ );
    //
    // path_to_tables_.swap( propagator.path_to_tables_ );

    // particle_->swap( *propagator.particle_ );
    // backup_particle_->swap( *propagator.backup_particle_ );
    //FirstOrderScattering
    // swap<ScatteringFirstOrder*> (scatteringFirstOrder_ ,propagator.scatteringFirstOrder_);
//    scatteringFirstOrderMoliere_->swap(*propagator.scatteringFirstOrderMoliere_);
    // swap(scattering_model_ , propagator.scattering_model_);

    // current_sector_->swap( *propagator.current_sector_ );
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//------------------------private member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// void Propagator::MoveParticle(double distance)
// {
//
//     double dist = particle_->GetPropagatedDistance();
//
//     Vector3D position = particle_->GetPosition();
//
//     dist   +=  distance;
//
//     position = position + distance*particle_->GetDirection();
//
//     particle_->SetPosition(position);
//
//     particle_->SetPropagatedDistance(dist);
//
// }


// void Propagator::InitDefaultCollection()
// {
//     Medium* med             = new Ice();
//     EnergyCutSettings* cuts = new EnergyCutSettings(500,0.05);
//     current_sector_ = new CollectionInterpolant(Ice(), Sphere(Vector3D(), 1e18, 0), EnergyCutSettings(500, 0.05));
//     current_sector_->SetGeometry(geom);
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// void Propagator::InitProcessCollections(ifstream &file)
// {
//
//     char buf[256];
//     string str_buf;
//     deque<string> *token = new deque<string>();
//     string taux;
//
//     unsigned int hirarchy;
//     bool found_hirarchy = false;
//
//     log_info("Reading sector informations");
//     double ecut_inside  = -1;
//     double ecut_infront = -1;
//     double ecut_behind  = -1;
//
//     double vcut_inside  = -1;
//     double vcut_infront = -1;
//     double vcut_behind  = -1;
//
//     bool cont_inside    =   false;
//     bool cont_infront   =   false;
//     bool cont_behind    =   false;
//
//     bool found_inside_cuts  =   false;
//     bool found_behind_cuts  =   false;
//     bool found_infront_cuts =   false;
//
//
//     while(file.good())
//     {
//         file.getline(buf,256);
//         str_buf =   string(buf);
//         token   =   SplitString(str_buf," \t");
//         if(!token->empty())
//         {
//             taux =   NextToken(token);
//             if(taux.at(0)!='#')
//             {
//                 break;
//             }
//         }
//     }
//     // Geometry *geometry  = new Geometry();
//     Geometry* geometry = InitGeometry(token,taux);
//     // detector_->SetHirarchy(0);
//     geometry->SetHirarchy(0);
//
//     while(file.good())
//     {
//         file.getline(buf,256);
//         str_buf =   string(buf);
//         token   =   SplitString(str_buf," \t");
//         if(token->empty())
//         {
//             continue;
//         }
//         taux =   NextToken(token);
//
//         if(taux.at(0)=='#')
//         {
//             continue;
//         }
//
//         if(ToLowerCase(taux).compare("hirarchy")==0)
//         {
//             found_hirarchy = true;
//             taux = NextToken(token);
//             try
//             {
//                 hirarchy  = boost::lexical_cast<unsigned int>(taux);
//             }
//             catch(boost::bad_lexical_cast&)
//             {
//                 log_warn("Chosen hirarchy %s but must be an unsigned int! Set to 0.", taux.c_str());
//                 hirarchy = 0;
//             }
//             geometry->SetHirarchy(hirarchy);
//             continue;
//         }
//         else if(ToLowerCase(taux).compare("inside")==0)
//         {
//             if(token->size()==3)
//             {
//                 taux    =   NextToken(token);
//                 try {
//                     ecut_inside  = boost::lexical_cast<double>(taux);
//                 }
//                 catch(boost::bad_lexical_cast&) {
//                     log_warn("ecut for inside of the detector is set to %s but must be a double! Set to 500.", taux.c_str());
//                     ecut_inside = 500;
//                 }
//
//                 taux    =   NextToken(token);
//                 try {
//                     vcut_inside  = boost::lexical_cast<double>(taux);
//                 }
//                 catch(boost::bad_lexical_cast&) {
//                     log_warn("vcut for inside of the detector is set to %s but must be a double! Set to -1.", taux.c_str());
//                     vcut_inside = -1;
//                 }
//
//                 taux    =   NextToken(token);
//                 try {
//                     cont_inside  = boost::lexical_cast<bool>(taux);
//                 }
//                 catch(boost::bad_lexical_cast&) {
//                     log_warn("cont for inside of the detector is set to %s but must be a bool! Set to false.", taux.c_str());
//                     cont_inside = false;
//                 }
//
//                 found_inside_cuts = true;
//
//             }
//             else
//             {
//                 log_warn("Expect 3 parameters afer keyword inside! Set inside cut settings to global cut settings");
//             }
//         }
//         else if(ToLowerCase(taux).compare("infront")==0)
//         {
//             if(token->size()==3)
//             {
//                 taux    =   NextToken(token);
//                 try {
//                     ecut_infront  = boost::lexical_cast<double>(taux);
//                 }
//                 catch(boost::bad_lexical_cast&) {
//                     log_warn("ecut for infront of the detector is set to %s but must be a double! Set to -1.", taux.c_str());
//                     ecut_infront = -1;
//                 }
//
//                 taux    =   NextToken(token);
//                 try {
//                     vcut_infront  = boost::lexical_cast<double>(taux);
//                 }
//                 catch(boost::bad_lexical_cast&) {
//                     log_warn("vcut for infront of the detector is set %s to but must be a double! Set to 0.001.", taux.c_str());
//                     vcut_infront = 0.001;
//                 }
//
//                 taux    =   NextToken(token);
//                 try {
//                     cont_infront  = boost::lexical_cast<bool>(taux);
//                 }
//                 catch(boost::bad_lexical_cast&) {
//                     log_warn("cont for infront of the detector is set to %s but must be a bool! Set to true.", taux.c_str());
//                     cont_infront = true;
//                 }
//
//                 found_infront_cuts = true;
//
//             }
//             else
//             {
//                 log_warn("Expect 3 parameters afer keyword infront! Set inside cut settings to global cut settings");
//             }
//         }
//         else if(ToLowerCase(taux).compare("behind")==0)
//         {
//             if(token->size()==3)
//             {
//                 taux    =   NextToken(token);
//                 try {
//                     ecut_behind  = boost::lexical_cast<double>(taux);
//                 }
//                 catch(boost::bad_lexical_cast&) {
//                     log_warn("ecut for behind of the detector is set to %s but must be a double! Set to -1.", taux.c_str());
//                     ecut_behind = -1;
//                 }
//
//                 taux    =   NextToken(token);
//                 try {
//                     vcut_behind  = boost::lexical_cast<double>(taux);
//                 }
//                 catch(boost::bad_lexical_cast&) {
//                     log_warn("vcut for behind of the detector is set to %s but must be a double! Set to -1.", taux.c_str());
//                     vcut_behind = -1;
//                 }
//
//                 taux    =   NextToken(token);
//                 try {
//                     cont_behind  = boost::lexical_cast<bool>(taux);
//                 }
//                 catch(boost::bad_lexical_cast&) {
//                     log_warn("cont for behind of the detector is set to %s but must be a bool! Set to false.", taux.c_str());
//                     cont_behind = false;
//                 }
//                 found_behind_cuts = true;
//
//             }
//             else
//             {
//                 log_warn("Expect 3 parameters afer keyword behind! Set inside cut settings to global cut settings");
//             }
//         }
//         else if(ToLowerCase(taux).compare("medium")==0)
//         {
//             string name     =   NextToken(token);
//
//             double density_correction;
//
//             taux    =   NextToken(token);
//             try {
//                 density_correction  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("density correction factor is set to %s but must be a double! Set to 1.",taux.c_str());
//                 density_correction = 1;
//             }
//
//             Medium* med = MediumFactory::Get()->CreateMedium(name);
//
//             if (particle_ == NULL)
//             {
//
//                 PROPOSALParticle *muminus    =   new PROPOSALParticle(MuMinusDef::Get());
//                 PROPOSALParticle *tauminus   =   new PROPOSALParticle(TauMinusDef::Get());
//                 PROPOSALParticle *eminus     =   new PROPOSALParticle(EMinusDef::Get());
//                 // PROPOSALParticle *muplus    =   new PROPOSALParticle(ParticleType::MuPlus);
//                 // PROPOSALParticle *tauplus   =   new PROPOSALParticle(ParticleType::TauPlus);
//                 // PROPOSALParticle *eplus     =   new PROPOSALParticle(ParticleType::EPlus);
//
//                 EnergyCutSettings *inside;
//                 EnergyCutSettings *infront;
//                 EnergyCutSettings *behind;
//
//                 ProcessCollection* muminus_inside;
//                 ProcessCollection* tauminus_inside;
//                 ProcessCollection* eminus_inside;
//                 // ProcessCollection* muplus_inside;
//                 // ProcessCollection* tauplus_inside;
//                 // ProcessCollection* eplus_inside;
//
//                 ProcessCollection* muminus_infront;
//                 ProcessCollection* tauminus_infront;
//                 ProcessCollection* eminus_infront;
//                 // ProcessCollection* muplus_infront;
//                 // ProcessCollection* tauplus_infront;
//                 // ProcessCollection* eplus_infront;
//
//                 ProcessCollection* muminus_behind;
//                 ProcessCollection* tauminus_behind;
//                 ProcessCollection* eminus_behind;
//                 // ProcessCollection* muplus_behind;
//                 // ProcessCollection* tauplus_behind;
//                 // ProcessCollection* eplus_behind;
//
//                 if(found_inside_cuts)
//                 {
//                     inside = new EnergyCutSettings(ecut_inside,vcut_inside);
//
//                     muminus_inside  = new ProcessCollection(new PROPOSALParticle(*muminus),med->clone(),new EnergyCutSettings(*inside));
//                     tauminus_inside = new ProcessCollection(new PROPOSALParticle(*tauminus),med->clone(),new EnergyCutSettings(*inside));
//                     eminus_inside   = new ProcessCollection(new PROPOSALParticle(*eminus),med->clone(),new EnergyCutSettings(*inside));
//                     // muplus_inside   = new ProcessCollection(new PROPOSALParticle(*muplus),new Medium(*med),new EnergyCutSettings(*inside));
//                     // tauplus_inside  = new ProcessCollection(new PROPOSALParticle(*tauplus),new Medium(*med),new EnergyCutSettings(*inside));
//                     // eplus_inside    = new ProcessCollection(new PROPOSALParticle(*eplus),new Medium(*med),new EnergyCutSettings(*inside));
//
//                     muminus_inside->SetEnableRandomization(cont_inside);
//                     muminus_inside->SetLocation(1);
//                     // muplus_inside->SetEnableRandomization(cont_inside);
//                     // muplus_inside->SetLocation(1);
//
//                     tauminus_inside->SetEnableRandomization(cont_inside);
//                     tauminus_inside->SetLocation(1);
//                     // tauplus_inside->SetEnableRandomization(cont_inside);
//                     // tauplus_inside->SetLocation(1);
//
//                     eminus_inside->SetEnableRandomization(cont_inside);
//                     eminus_inside->SetLocation(1);
//                     // eplus_inside->SetEnableRandomization(cont_inside);
//                     // eplus_inside->SetLocation(1);
//                 }
//                 else
//                 {
//                     inside = new EnergyCutSettings(global_ecut_inside_,global_vcut_inside_);
//
//                     muminus_inside  = new ProcessCollection(new PROPOSALParticle(*muminus),med->clone(),new EnergyCutSettings(*inside));
//                     tauminus_inside = new ProcessCollection(new PROPOSALParticle(*tauminus),med->clone(),new EnergyCutSettings(*inside));
//                     eminus_inside   = new ProcessCollection(new PROPOSALParticle(*eminus),med->clone(),new EnergyCutSettings(*inside));
//                     // muplus_inside   = new ProcessCollection(new PROPOSALParticle(*muplus),new Medium(*med),new EnergyCutSettings(*inside));
//                     // tauplus_inside  = new ProcessCollection(new PROPOSALParticle(*tauplus),new Medium(*med),new EnergyCutSettings(*inside));
//                     // eplus_inside    = new ProcessCollection(new PROPOSALParticle(*eplus),new Medium(*med),new EnergyCutSettings(*inside));
//
//                     muminus_inside->SetEnableRandomization(global_cont_inside_);
//                     muminus_inside->SetLocation(1);
//                     // muplus_inside->SetEnableRandomization(global_cont_inside_);
//                     // muplus_inside->SetLocation(1);
//
//                     tauminus_inside->SetEnableRandomization(global_cont_inside_);
//                     tauminus_inside->SetLocation(1);
//                     // tauplus_inside->SetEnableRandomization(global_cont_inside_);
//                     // tauplus_inside->SetLocation(1);
//
//                     eminus_inside->SetEnableRandomization(global_cont_inside_);
//                     eminus_inside->SetLocation(1);
//                     // eplus_inside->SetEnableRandomization(global_cont_inside_);
//                     // eplus_inside->SetLocation(1);
//                 }
//
//                 if(found_infront_cuts)
//                 {
//                     infront = new EnergyCutSettings(ecut_infront,vcut_infront);
//
//                     muminus_infront  = new ProcessCollection(new PROPOSALParticle(*muminus),med->clone(),new EnergyCutSettings(*infront));
//                     tauminus_infront = new ProcessCollection(new PROPOSALParticle(*tauminus),med->clone(),new EnergyCutSettings(*infront));
//                     eminus_infront   = new ProcessCollection(new PROPOSALParticle(*eminus),med->clone(),new EnergyCutSettings(*infront));
//                     // muplus_infront   = new ProcessCollection(new PROPOSALParticle(*muplus),new Medium(*med),new EnergyCutSettings(*infront));
//                     // tauplus_infront  = new ProcessCollection(new PROPOSALParticle(*tauplus),new Medium(*med),new EnergyCutSettings(*infront));
//                     // eplus_infront    = new ProcessCollection(new PROPOSALParticle(*eplus),new Medium(*med),new EnergyCutSettings(*infront));
//
//                     muminus_infront->SetEnableRandomization(cont_infront);
//                     muminus_infront->SetLocation(0);
//                     // muplus_infront->SetEnableRandomization(cont_infront);
//                     // muplus_infront->SetLocation(0);
//
//                     tauminus_infront->SetEnableRandomization(cont_infront);
//                     tauminus_infront->SetLocation(0);
//                     // tauplus_infront->SetEnableRandomization(cont_infront);
//                     // tauplus_infront->SetLocation(0);
//
//                     eminus_infront->SetEnableRandomization(cont_infront);
//                     eminus_infront->SetLocation(0);
//                     // eplus_infront->SetEnableRandomization(cont_infront);
//                     // eplus_infront->SetLocation(0);
//                 }
//                 else
//                 {
//                     infront = new EnergyCutSettings(global_ecut_infront_,global_vcut_infront_);
//
//                     muminus_infront  = new ProcessCollection(new PROPOSALParticle(*muminus),med->clone(),new EnergyCutSettings(*infront));
//                     tauminus_infront = new ProcessCollection(new PROPOSALParticle(*tauminus),med->clone(),new EnergyCutSettings(*infront));
//                     eminus_infront   = new ProcessCollection(new PROPOSALParticle(*eminus),med->clone(),new EnergyCutSettings(*infront));
//                     // muplus_infront   = new ProcessCollection(new PROPOSALParticle(*muplus),new Medium(*med),new EnergyCutSettings(*infront));
//                     // tauplus_infront  = new ProcessCollection(new PROPOSALParticle(*tauplus),new Medium(*med),new EnergyCutSettings(*infront));
//                     // eplus_infront    = new ProcessCollection(new PROPOSALParticle(*eplus),new Medium(*med),new EnergyCutSettings(*infront));
//
//                     muminus_infront->SetEnableRandomization(global_cont_infront_);
//                     muminus_infront->SetLocation(0);
//                     // muplus_infront->SetEnableRandomization(global_cont_infront_);
//                     // muplus_infront->SetLocation(0);
//
//                     tauminus_infront->SetEnableRandomization(global_cont_infront_);
//                     tauminus_infront->SetLocation(0);
//                     // tauplus_infront->SetEnableRandomization(global_cont_infront_);
//                     // tauplus_infront->SetLocation(0);
//
//                     eminus_infront->SetEnableRandomization(global_cont_infront_);
//                     eminus_infront->SetLocation(0);
//                     // eplus_infront->SetEnableRandomization(global_cont_infront_);
//                     // eplus_infront->SetLocation(0);
//                 }
//
//                 if(found_behind_cuts)
//                 {
//                     behind = new EnergyCutSettings(ecut_behind,vcut_behind);
//
//                     muminus_behind  = new ProcessCollection(new PROPOSALParticle(*muminus),med->clone(),new EnergyCutSettings(*behind));
//                     tauminus_behind = new ProcessCollection(new PROPOSALParticle(*tauminus),med->clone(),new EnergyCutSettings(*behind));
//                     eminus_behind   = new ProcessCollection(new PROPOSALParticle(*eminus),med->clone(),new EnergyCutSettings(*behind));
//                     // muplus_behind   = new ProcessCollection(new PROPOSALParticle(*muplus),new Medium(*med),new EnergyCutSettings(*behind));
//                     // tauplus_behind  = new ProcessCollection(new PROPOSALParticle(*tauplus),new Medium(*med),new EnergyCutSettings(*behind));
//                     // eplus_behind    = new ProcessCollection(new PROPOSALParticle(*eplus),new Medium(*med),new EnergyCutSettings(*behind));
//
//                     muminus_behind->SetEnableRandomization(cont_behind);
//                     muminus_behind->SetLocation(2);
//                     // muplus_behind->SetEnableRandomization(cont_behind);
//                     // muplus_behind->SetLocation(2);
//
//                     tauminus_behind->SetEnableRandomization(cont_behind);
//                     tauminus_behind->SetLocation(2);
//                     // tauplus_behind->SetEnableRandomization(cont_behind);
//                     // tauplus_behind->SetLocation(2);
//
//                     eminus_behind->SetEnableRandomization(cont_behind);
//                     eminus_behind->SetLocation(2);
//                     // eplus_behind->SetEnableRandomization(cont_behind);
//                     // eplus_behind->SetLocation(2);
//                 }
//                 else
//                 {
//                     behind = new EnergyCutSettings(global_ecut_behind_,global_vcut_behind_);
//
//                     muminus_behind  = new ProcessCollection(new PROPOSALParticle(*muminus),med->clone(),new EnergyCutSettings(*behind));
//                     tauminus_behind = new ProcessCollection(new PROPOSALParticle(*tauminus),med->clone(),new EnergyCutSettings(*behind));
//                     eminus_behind   = new ProcessCollection(new PROPOSALParticle(*eminus),med->clone(),new EnergyCutSettings(*behind));
//                     // muplus_behind   = new ProcessCollection(new PROPOSALParticle(*muplus),new Medium(*med),new EnergyCutSettings(*behind));
//                     // tauplus_behind  = new ProcessCollection(new PROPOSALParticle(*tauplus),new Medium(*med),new EnergyCutSettings(*behind));
//                     // eplus_behind    = new ProcessCollection(new PROPOSALParticle(*eplus),new Medium(*med),new EnergyCutSettings(*behind));
//
//                     muminus_behind->SetEnableRandomization(global_cont_behind_);
//                     muminus_behind->SetLocation(2);
//                     // muplus_behind->SetEnableRandomization(global_cont_behind_);
//                     // muplus_behind->SetLocation(2);
//
//                     tauminus_behind->SetEnableRandomization(global_cont_behind_);
//                     tauminus_behind->SetLocation(2);
//                     // tauplus_behind->SetEnableRandomization(global_cont_behind_);
//                     // tauplus_behind->SetLocation(2);
//
//                     eminus_behind->SetEnableRandomization(global_cont_behind_);
//                     eminus_behind->SetLocation(2);
//                     // eplus_behind->SetEnableRandomization(global_cont_behind_);
//                     // eplus_behind->SetLocation(2);
//                 }
//
//
//                 int former_size =sectors_.size();
//
//                 sectors_.push_back( muminus_infront );
//                 sectors_.push_back( muminus_inside );
//                 sectors_.push_back( muminus_behind );
//                 // sectors_.push_back( muplus_infront );
//                 // sectors_.push_back( muplus_inside );
//                 // sectors_.push_back( muplus_behind );
//
//                 sectors_.push_back( tauminus_infront );
//                 sectors_.push_back( tauminus_inside );
//                 sectors_.push_back( tauminus_behind );
//                 // sectors_.push_back( tauplus_infront );
//                 // sectors_.push_back( tauplus_inside );
//                 // sectors_.push_back( tauplus_behind );
//
//                 sectors_.push_back( eminus_infront );
//                 sectors_.push_back( eminus_inside );
//                 sectors_.push_back( eminus_behind );
//                 // sectors_.push_back( eplus_infront );
//                 // sectors_.push_back( eplus_inside );
//                 // sectors_.push_back( eplus_behind );
//
//                 for(unsigned int i = former_size ;i<sectors_.size(); i++)
//                 {
//                     sectors_.at(i)->SetGeometry(geometry);
//                     sectors_.at(i)->SetDensityCorrection(density_correction);
//                 }
//
//                 delete med;
//                 delete muminus;
//                 delete tauminus;
//                 delete eminus;
//                 // delete muplus;
//                 // delete tauplus;
//                 // delete eplus;
//                 delete inside;
//                 delete infront;
//                 delete behind;
//
//                 break;
//             }
//             else
//             {
//                 log_info("Got particle from constructor:%s with mass=%f", particle_->GetName().c_str(), particle_->GetMass());
//                 EnergyCutSettings* inside;
//                 EnergyCutSettings* infront;
//                 EnergyCutSettings* behind;
//
//                 ProcessCollection* particle_inside;
//                 ProcessCollection* particle_infront;
//                 ProcessCollection* particle_behind;
//
//                 if(found_inside_cuts)
//                 {
//                     inside = new EnergyCutSettings(ecut_inside,vcut_inside);
//                     particle_inside  = new ProcessCollection(new PROPOSALParticle(*particle_),med->clone(),new EnergyCutSettings(*inside));
//                     particle_inside->SetEnableRandomization(cont_inside);
//                     particle_inside->SetLocation(1);
//                 }
//                 else
//                 {
//                     inside = new EnergyCutSettings(global_ecut_inside_,global_vcut_inside_);
//                     particle_inside  = new ProcessCollection(new PROPOSALParticle(*particle_),med->clone(),new EnergyCutSettings(*inside));
//                     particle_inside->SetEnableRandomization(global_cont_inside_);
//                     particle_inside->SetLocation(1);
//                 }
//
//                 if (found_infront_cuts)
//                 {
//                     infront = new EnergyCutSettings(ecut_infront,vcut_infront);
//                     particle_infront  = new ProcessCollection(new PROPOSALParticle(*particle_),med->clone(),new EnergyCutSettings(*infront));
//                     particle_infront->SetEnableRandomization(cont_infront);
//                     particle_infront->SetLocation(0);
//                 }
//                 else
//                 {
//                     infront = new EnergyCutSettings(global_ecut_infront_,global_vcut_infront_);
//                     particle_infront  = new ProcessCollection(new PROPOSALParticle(*particle_),med->clone(),new EnergyCutSettings(*infront));
//                     particle_infront->SetEnableRandomization(global_cont_infront_);
//                     particle_infront->SetLocation(0);
//                 }
//
//                 if(found_behind_cuts)
//                 {
//                     behind = new EnergyCutSettings(ecut_behind,vcut_behind);
//                     particle_behind  = new ProcessCollection(new PROPOSALParticle(*particle_),med->clone(),new EnergyCutSettings(*behind));
//                     particle_behind->SetEnableRandomization(cont_behind);
//                     particle_behind->SetLocation(2);
//                 }
//                 else
//                 {
//                     behind = new EnergyCutSettings(global_ecut_behind_,global_vcut_behind_);
//                     particle_behind  = new ProcessCollection(new PROPOSALParticle(*particle_),med->clone(),new EnergyCutSettings(*behind));
//                     particle_behind->SetEnableRandomization(global_cont_behind_);
//                     particle_behind->SetLocation(2);
//                 }
//
//                 int former_size =sectors_.size();
//                 sectors_.push_back( particle_infront );
//                 sectors_.push_back( particle_inside );
//                 sectors_.push_back( particle_behind );
//
//                 for(unsigned int i = former_size ;i<sectors_.size(); i++)
//                 {
//                     sectors_.at(i)->SetGeometry(geometry);
//                     sectors_.at(i)->SetDensityCorrection(density_correction);
//                 }
//
//                 delete med;
//                 delete inside;
//                 delete infront;
//                 delete behind;
//                 break;
//             }
//         }
//         else
//         {
//             log_fatal("Last line in a sector segment must start with ’medium’!");
//             exit(1);
//         }
//     }
//
//     if(!found_hirarchy)log_info("Hirarchy for geometry not set!");
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//



int Propagator::GetSeed() const
{
    return seed_;
}

void Propagator::SetSeed(int seed)
{
    seed_ = seed;
}

// ParametrizationType::Enum Propagator::GetBrems() const
// {
//     return brems_;
// }
//
// void Propagator::SetBrems(ParametrizationType::Enum brems)
// {
//     brems_ = brems;
// }
//
// ParametrizationType::Enum Propagator::GetPhoto() const
// {
//     return photo_;
// }
//
// void Propagator::SetPhoto(ParametrizationType::Enum photo)
// {
//     photo_ = photo;
// }

// std::string Propagator::GetPath_to_tables() const
// {
//     return path_to_tables_;
// }
//
// void Propagator::SetPath_to_tables(const std::string &path_to_tables)
// {
//     path_to_tables_ = path_to_tables;
// }

Geometry *Propagator::GetDetector() const
{
    return detector_;
}

PROPOSALParticle& Propagator::GetParticle()
{
    return particle_;
}

// void Propagator::SetDetector(Geometry *detector)
// {
//     detector_ = detector;
// }

// bool Propagator::GetStopping_decay() const
// {
//     return stopping_decay_;
// }
//
// void Propagator::SetStopping_decay(bool stopping_decay)
// {
//     stopping_decay_ = stopping_decay;
// }

// void Propagator::RestoreBackup_particle()
// {
//     particle_ = new PROPOSALParticle(*backup_particle_);
// }

// void Propagator::ResetParticle()
// {
//     // particle_ = new PROPOSALParticle(*backup_particle_);
//     *particle_ = *backup_particle_;
// }

// Geometry* Propagator::InitGeometry(std::deque<std::string>* token, string first_token)
// {
//     string taux = first_token;
//
//     double x_coord = 0.;
//     double y_coord = 0.;
//     double z_coord = 0.;
//
//     if(ToLowerCase(taux).compare("cylinder")==0)
//     {
//         double radius;
//         double inner_radius =   0;
//         double height;
//
//         // Only radius and height are specified
//         if(token->size()==2)
//         {
//             taux    =   NextToken(token);
//             try {
//                 radius  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_fatal("radius must be a double! Exit");
//                 exit(1);
//             }
//
//             taux    =   NextToken(token);
//             try {
//                 height  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_fatal("height must be a double! Exit");
//                 exit(1);
//             }
//         }
//         // Only radius, inner_radius and height are specified
//         else if(token->size()==3)
//         {
//             taux    =   NextToken(token);
//             try {
//                 radius  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_fatal("adius must be a double! Exit");
//                 exit(1);
//             }
//
//             taux    =   NextToken(token);
//             try {
//                 inner_radius  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("inner_radius is set to %s but must be a double! Set to 0",taux.c_str());
//                 inner_radius    =   0;
//             }
//
//             taux    =   NextToken(token);
//             try {
//                 height  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_fatal("Error: height must be a double! Exit");
//                 exit(1);
//             }
//         }
//         // position vector, radius, inner_radius and height are specified
//         else if(token->size()==6)
//         {
//             taux    =   NextToken(token);
//             try {
//                 x_coord = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("x position is set to %s but must be a double! Set to 0",taux.c_str());
//                 x_coord = 0;
//             }
//
//             taux    =   NextToken(token);
//             try {
//                 y_coord = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("y position is set to %s but must be a double! Set to 0", taux.c_str());
//                 y_coord = 0;
//             }
//
//             taux    =   NextToken(token);
//             try {
//                 z_coord = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("z position is set to %s but must be a double! Set to 0",taux.c_str());
//                 z_coord = 0;
//             }
//
//             taux    =   NextToken(token);
//             try {
//                 radius  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_fatal("radius must be a double! Exit");
//                 exit(1);
//             }
//
//             taux    =   NextToken(token);
//             try {
//                 inner_radius  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("inner_radius is set to %s but must be a double! Set to 0", taux.c_str());
//                 inner_radius    =   0;
//             }
//
//             taux    =   NextToken(token);
//             try {
//                 height  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_fatal("Error: height must be a double! Exit");
//                 exit(1);
//             }
//         }
//         else
//         {
//             log_fatal("Number of values after 'cylinder' must be 2,3 or 6. Exit!");
//             exit(1);
//         }
//
//         return new Cylinder(Vector3D(x_coord, y_coord, z_coord),radius,inner_radius,height);
//         // geometry->InitCylinder(origin_x,origin_y,origin_z,radius,inner_radius,height);
//
//     }
//     else if(ToLowerCase(taux).compare("sphere")==0)
//     {
//         double radius;
//         double inner_radius =   0;
//
//         // Only radius is specified
//         if(token->size()==1)
//         {
//             taux    =   NextToken(token);
//             try {
//                 radius  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_fatal("radius must be a double! Exit");
//                 exit(1);
//             }
//         }
//         // Only radius and inner_radius are specified
//         else if(token->size()==2)
//         {
//             taux    =   NextToken(token);
//             try {
//                 radius  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_fatal("radius must be a double! Exit");
//                 exit(1);
//             }
//
//             taux    =   NextToken(token);
//             try {
//                 inner_radius  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("inner_radius is set to %s but must be a double! Set to 0", taux.c_str());
//                 inner_radius    =   0;
//             }
//
//         }
//         // position vector, radius, inner_radius and height are specified
//         else if(token->size()==5)
//         {
//             taux    =   NextToken(token);
//             try {
//                 x_coord = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("x position is set to %s but must be a double! Set to 0", taux.c_str());
//                 x_coord = 0;
//             }
//
//             taux    =   NextToken(token);
//             try {
//                 y_coord = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("y position is set to %s but must be a double! Set to 0", taux.c_str());
//                 y_coord = 0;
//             }
//
//             taux    =   NextToken(token);
//             try {
//                 z_coord = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("z position is set to %s but must be a double! Set to 0", taux.c_str());
//                 z_coord = 0;
//             }
//
//             taux    =   NextToken(token);
//             try {
//                 radius  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_fatal("Error: radius must be a double! Exit");
//                 exit(1);
//             }
//
//             taux    =   NextToken(token);
//             try {
//                 inner_radius  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("inner_radius is set to %s but must be a double! Set to 0", taux.c_str());
//                 inner_radius    =   0;
//             }
//
//         }
//         else
//         {
//             log_fatal("Number of values after 'sphere' must be 1,2 or 5. Exit!");
//             exit(1);
//         }
//
//         return new Sphere(Vector3D(x_coord, y_coord, z_coord),radius,inner_radius);
//         // geometry->InitSphere(origin_x,origin_y,origin_z,radius,inner_radius);
//
//     }
//     else if(ToLowerCase(taux).compare("box")==0)
//     {
//         double width_x;
//         double width_y;
//         double width_z;
//
//         // Only width_x,width_y and width_z are specified
//         if(token->size()==3)
//         {
//             taux    =   NextToken(token);
//             try {
//                 width_x  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_fatal("width_x must be a double! Exit");
//                 exit(1);
//             }
//
//             taux    =   NextToken(token);
//             try {
//                 width_y  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_fatal("width_y must be a double! Exit");
//                 exit(1);
//             }
//
//             taux    =   NextToken(token);
//             try {
//                 width_z  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_fatal("width_z must be a double! Exit");
//                 exit(1);
//             }
//         }
//         // position vector, width_x,width_y and width_z are specified
//         else if(token->size()==6)
//         {
//             taux    =   NextToken(token);
//             try {
//                 x_coord = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("x position is set to %s but must be a double! Set to 0", taux.c_str());
//                 x_coord = 0;
//             }
//
//             taux    =   NextToken(token);
//             try {
//                 y_coord = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("y position is set to %s but must be a double! Set to 0", taux.c_str());
//                 y_coord = 0;
//             }
//
//             taux    =   NextToken(token);
//             try {
//                 z_coord = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_warn("z position is set to %s but must be a double! Set to 0", taux.c_str());
//                 z_coord = 0;
//             }
//
//             taux    =   NextToken(token);
//             try {
//                 width_x  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_fatal("width_x must be a double! Exit");
//                 exit(1);
//             }
//
//             taux    =   NextToken(token);
//             try {
//                 width_y  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_fatal("width_y must be a double! Exit");
//                 exit(1);
//             }
//
//             taux    =   NextToken(token);
//             try {
//                 width_z  = boost::lexical_cast<double>(taux);
//             }
//             catch(boost::bad_lexical_cast&) {
//                 log_fatal("width_z must be a double! Exit");
//                 exit(1);
//             }
//
//         }
//         else
//         {
//             log_fatal("Number of values after 'box' must be 3 or 6. Exit!");
//             exit(1);
//         }
//
//         return new Box(Vector3D(x_coord, y_coord, z_coord),width_x,width_y,width_z);
//         // geometry->InitBox(origin_x,origin_y,origin_z,width_x,width_y,width_z);
//
//     }
//     else
//     {
//         log_fatal("Unrecognized geometry: %s Must be cylinder, sphere or box! Exit",taux.c_str());
//         exit(1);
//     }
//
// }
