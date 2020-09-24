/*
 * Propagator.cxx
 * Created on: 23.04.2013
 * Author: koehne
 */

// #include <cmath>

#include <fstream>
#include <memory>

#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"

#include "PROPOSAL/geometry/Box.h"
#include "PROPOSAL/geometry/Cylinder.h"
#include "PROPOSAL/geometry/GeometryFactory.h"
#include "PROPOSAL/geometry/Sphere.h"

#include "PROPOSAL/particle/Particle.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
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
    const std::vector<Sector*>& sectors, std::shared_ptr<const Geometry> geometry) try
    : current_sector_(NULL),
      particle_def_(sectors.at(0)->GetParticleDef()),
      detector_(geometry)
{
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
    std::shared_ptr<const Geometry> geometry)
    : particle_def_(particle_def)
    , detector_(geometry)
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
    std::shared_ptr<const Geometry> geometry, const InterpolationDef& interpolation_def)
    : particle_def_(particle_def)
    , detector_(geometry)
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
    , detector_(propagator.detector_)
{
    for (unsigned int i = 0; i < propagator.sectors_.size(); ++i) {
        sectors_[i] = new Sector(*propagator.sectors_[i]);

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
{
    int global_seed = global_seed_;
    bool do_interpolation = do_interpolation_;
    bool uniform = uniform_;

    InterpolationDef interpolation_def;
    std::unique_ptr<Sector::Definition> sec_def_global(new Sector::Definition());

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

    nlohmann::json json_global;
    if(json_config.contains("global")){
        json_global = json_config["global"];

        global_seed = json_global.value("seed", global_seed_);
        uniform = json_global.value("uniform", uniform_);

        if (json_global.contains("interpolation")){
            nlohmann::json json_interpol = json_global["interpolation"];
            interpolation_def = InterpolationDef(json_interpol);
            do_interpolation = json_interpol.value("do_interpolation", do_interpolation_);
        }

        sec_def_global.reset(new Sector::Definition(json_global));
    }

    RandomGenerator::Get().SetSeed(global_seed);
    particle_def_.decay_table.SetUniformSampling(uniform);

    if (json_config.find("detector") != json_config.end()) {
        std::string shape = json_config["detector"]["shape"];
        if (shape == "sphere") {
            detector_ = std::make_shared<const Sphere>(json_config["detector"]);
        } else if (shape == "box") {
            detector_ = std::make_shared<const Box>(json_config["detector"]);
        } else if (shape == "cylinder") {
            detector_ = std::make_shared<const Cylinder>(json_config["detector"]);
        } else {
            throw std::invalid_argument("You need to specify a detector for each sector");
        }
    }

    std::shared_ptr<const Medium> med;
    std::shared_ptr<const Geometry> geo;
    std::array<std::pair<std::string, Sector::ParticleLocation::Enum>, 3> cuts {
        std::make_pair("cuts_infront", Sector::ParticleLocation::InfrontDetector),
        std::make_pair("cuts_inside", Sector::ParticleLocation::InsideDetector),
        std::make_pair("cuts_behind", Sector::ParticleLocation::BehindDetector)
    };
    for (const auto& cut : cuts) {
        if(json_global.contains(cut.first)) {
            sec_def_global->cut_settings = EnergyCutSettings(json_global.at(cut.first));
            sec_def_global->location = cut.second;
        }

        if (json_config.contains("sectors")) {
            assert(json_config["sectors"].is_array());
            for (const auto& json_sector : json_config.at("sectors")) {
                double density_correction = json_sector.value("density_correction", 1.0);
                std::string medium_name = json_sector.value("medium","water");
                med = CreateMedium(medium_name, density_correction);
                if (json_sector.contains("geometry"))
                {
                    std::string shape = json_sector["geometry"]["shape"];
                    if (shape == "sphere") {
                        geo = std::make_shared<const Sphere>(json_sector["geometry"]);
                    } else if (shape == "box") {
                        geo = std::make_shared<const Box>(json_sector["geometry"]);
                    } else if (shape == "cylinder") {
                        geo = std::make_shared<const Cylinder>(json_sector["geometry"]);
                    } else {
                        throw std::invalid_argument("You need to specify a detector for each sector");
                    }
                }

                    Sector::Definition sec_def = *sec_def_global;
                    sec_def.SetMedium(med);
                    sec_def.SetGeometry(geo);

                    if(json_sector.contains(cut.first))
                    {
                        nlohmann::json local_cuts = json_sector[cut.first];
                        sec_def.cut_settings = EnergyCutSettings(local_cuts);
                        if(local_cuts.contains("cont_rand")) {
                            sec_def.do_continuous_randomization = local_cuts.value("cont_rand", true);
                        }
                    }

                    if (do_interpolation) {
                        sectors_.push_back(new Sector(particle_def_, sec_def, interpolation_def));
                    } else {
                        sectors_.push_back(new Sector(particle_def_, sec_def));
                    }
                }
            }
    }
}

Propagator::~Propagator()
{
    for (auto sector : sectors_) {
        delete sector;
    }

    sectors_.clear();
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
    const DynamicData& initial_condition, double max_distance, double minimal_energy)
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

    Vector3D position(initial_condition.GetPosition());
    Vector3D direction(initial_condition.GetDirection());

    /* bool starts_in_detector = detector_->IsInside( initial_condition.GetPosition(), initial_condition.GetDirection()); */
    bool starts_in_detector = detector_->IsInside( position, direction);
    if (starts_in_detector) {
        secondaries_.SetEntryPoint(initial_condition);
        distance_to_closest_approach = detector_->DistanceToClosestApproach(
            initial_condition.GetPosition(), initial_condition.GetDirection());
        if (distance_to_closest_approach < 0) {
            secondaries_.SetClosestApproachPoint(initial_condition);
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
                        secondaries_.SetClosestApproachPoint(*p_condition);
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
            secondaries_.SetClosestApproachPoint(*p_condition);

            propagationstep_till_closest_approach = false;
        }

        if (std::abs(max_distance - p_condition->GetPropagatedDistance()) < PARTICLE_POSITION_RESOLUTION
            || p_condition->GetEnergy() <= minimal_energy
            || p_condition->GetType() == static_cast<int>(InteractionType::Decay))
            break;
    }
    if (detector_->IsInside(
            p_condition->GetPosition(), p_condition->GetDirection())) {
        secondaries_.SetExitPoint(*p_condition);
    }

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

    // Get Location of the detector (Inside/Infront/Behind)
    Geometry::ParticleLocation::Enum detector_location
        = detector_->GetLocation(particle_position, particle_direction);
    for (unsigned int i = 0; i < sectors_.size(); ++i) {
        if (sectors_[i]->GetSectorDef().GetGeometry()->IsInside(particle_position, particle_direction)) {
            if (static_cast<int>(sectors_[i]->GetLocation()) == static_cast<int>(detector_location))
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
        if (current_sector_->GetSectorDef().GetGeometry()->GetHierarchy()
            == sectors_[*iter]->GetSectorDef().GetGeometry()->GetHierarchy()) {
            // Current Density is smaller -> Set the new sector!
            //
            if (current_sector_->GetSectorDef().GetMedium()->GetCorrectedMassDensity(
                    particle_position)
                < sectors_[*iter]->GetSectorDef().GetMedium()->GetCorrectedMassDensity(
                      particle_position))
                current_sector_ = sectors_[*iter];
        }

        // Current Hierarchy is smaller -> Set the new sector!
        //
        if (current_sector_->GetSectorDef().GetGeometry()->GetHierarchy()
            < sectors_[*iter]->GetSectorDef().GetGeometry()->GetHierarchy())
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
        = current_sector_->GetSectorDef().GetGeometry()
              ->DistanceToBorder(particle_position, particle_direction)
              .first;
    double tmp_distance_to_border;

    Geometry::ParticleLocation::Enum detector_location
        = detector_->GetLocation(particle_position, particle_direction);

    for (std::vector<Sector*>::iterator iter = sectors_.begin();
         iter != sectors_.end(); ++iter) {

        if (static_cast<int>((*iter)->GetLocation())
            == static_cast<int>(detector_location)) {
            if ((*iter)->GetSectorDef().GetGeometry()->GetHierarchy()
                >= current_sector_->GetSectorDef().GetGeometry()->GetHierarchy()) {
                tmp_distance_to_border
                    = (*iter)
                          ->GetSectorDef().GetGeometry()
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
