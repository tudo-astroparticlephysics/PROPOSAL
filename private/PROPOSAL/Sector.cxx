
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/Secondaries.h"

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/decay/DecayChannel.h"

#include "PROPOSAL/medium/density_distr/density_distr.h"

#include "PROPOSAL/particle/Particle.h"

#include "PROPOSAL/geometry/Geometry.h"
#include "PROPOSAL/geometry/GeometryFactory.h"
#include "PROPOSAL/geometry/Sphere.h"

#include "PROPOSAL/Sector.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"

#include "PROPOSAL/methods.h"
#include "PROPOSAL/propagation_utility/ContinuousRandomizer.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"

#include <algorithm>
#include <array>
#include <memory>
#include <utility>

using namespace PROPOSAL;

namespace PROPOSAL {

std::ostream& operator<<(std::ostream& os, Sector::Definition const& sec_definition)
{
    std::stringstream ss;
    ss << " Sector Definition (" << &sec_definition << ") ";
    os << Helper::Centered(60, ss.str()) << '\n';

    os << "Do Stochastic Loss Weighting: " << sec_definition.do_stochastic_loss_weighting << std::endl;
    os << "Stochastic Loss Weighting: " << sec_definition.stochastic_loss_weighting << std::endl;
    os << "Dp Stopping Decay: " << sec_definition.stopping_decay << std::endl;
    os << "Dp Continuous Randomization: " << sec_definition.do_continuous_randomization << std::endl;
    os << "Dp Continuous Energy Loss output: " << sec_definition.do_continuous_energy_loss_output << std::endl;
    os << "Dp Exact Time calculation: " << sec_definition.do_exact_time_calculation << std::endl;
    os << "Only store loss inside the detector volume: " << sec_definition.only_loss_inside_detector << std::endl;
    os << "Scattering Model: " << sec_definition.scattering_model << std::endl;
    os << "Particle Location: " << sec_definition.location << std::endl;
    os << "Propagation Utility Definition:\n" << sec_definition.utility_def << std::endl;
    os << "Energy Cut Settings:\n" << sec_definition.cut_settings << std::endl;
    os << "Medium:\n" << *sec_definition.medium_ << std::endl;
    os << "Geometry:\n" << *sec_definition.geometry_ << std::endl;

    os << Helper::Centered(60, "");
    return os;
}

} // namespace PROPOSAL

/******************************************************************************
 *                                 Sector                                 *
 ******************************************************************************/

Sector::Definition::Definition()
    : do_stochastic_loss_weighting(false)
    , stochastic_loss_weighting(0)
    , stopping_decay(true)
    , do_continuous_randomization(true)
    , do_continuous_energy_loss_output(false)
    , do_exact_time_calculation(true)
    , only_loss_inside_detector(false)
    , scattering_model(ScatteringFactory::HighlandIntegral)
    , location(Sector::ParticleLocation::InsideDetector)
    , utility_def()
    , cut_settings()
    , medium_(std::make_shared<const Ice>())
    , geometry_(std::make_shared<const Sphere>(Vector3D(), 1.0e20, 0.0))
{
}

Sector::Definition::Definition(const Definition& def)
    : do_stochastic_loss_weighting(def.do_stochastic_loss_weighting)
    , stochastic_loss_weighting(def.stochastic_loss_weighting)
    , stopping_decay(def.stopping_decay)
    , do_continuous_randomization(def.do_continuous_randomization)
    , do_continuous_energy_loss_output(def.do_continuous_energy_loss_output)
    , do_exact_time_calculation(def.do_exact_time_calculation)
    , only_loss_inside_detector(def.only_loss_inside_detector)
    , scattering_model(def.scattering_model)
    , location(def.location)
    , utility_def(def.utility_def)
    , cut_settings(def.cut_settings)
    , medium_(def.medium_)
    , geometry_(def.geometry_)
{
}

Sector::Definition::Definition(const nlohmann::json& config)
    : utility_def(config)
{
    if(not config.is_object()) throw std::invalid_argument("No json object found.");

    std::array<std::pair<std::string, Sector::ParticleLocation::Enum>, 3> cuts = {
        std::make_pair("cuts_infront", Sector::ParticleLocation::InfrontDetector),
        std::make_pair("cuts_inside",Sector::ParticleLocation::InsideDetector),
        std::make_pair("cuts_behind",Sector::ParticleLocation::BehindDetector),
    };

    for (const auto& cut : cuts) {
        if(config.contains(cut.first)) {
            cut_settings = EnergyCutSettings(config.at(cut.first));
            location = cut.second;
            break;
        }
    }

    std::string scattering_model_str = config.value("scattering", "highlandintegral");
    scattering_model = ScatteringFactory::Get().GetEnumFromString(scattering_model_str);

    if(!config.contains("geometry")) {
        throw std::invalid_argument("A geometry is needed for each sector.");
    }
    if(!config.at("geometry").contains("shape")) {
        throw std::invalid_argument("A shape must be specified.");
    }

    nlohmann::json config_geometry = config.at("geometry");
    if (config_geometry.at("shape") == "sphere") {
        geometry_ = std::make_shared<const Sphere>(config_geometry);
    } else if (config_geometry.at("shape") == "box") {
        geometry_ = std::make_shared<const Box>(config_geometry);
    } else if (config_geometry.at("shape") == "cylinder") {
        geometry_ = std::make_shared<const Cylinder>(config_geometry);
    } else {
        throw std::invalid_argument("Unkown shape.");
    }

    if(!config.contains("medium"))
    {
        throw std::invalid_argument("A medium is needed for each sector.");
    }

    std::string medium_name;
    config.at("medium").get_to(medium_name);
    double density_correction = config.value("density_correction", 1.);
    medium_ = CreateMedium(medium_name, density_correction);


    do_stochastic_loss_weighting = config.value("do_stochastic_loss_weighting", false);
    stochastic_loss_weighting = config.value("stochastic_loss_weighting", 0);
    stopping_decay = config.value("stopping_decay", true);
    do_continuous_randomization = config.value("cont_rand", true);
    do_continuous_energy_loss_output = config.value("do_continuous_energy_loss_output", false);
    do_exact_time_calculation = config.value("exact_time", true);
    only_loss_inside_detector = config.value("only_loss_inside_detector", false);
}

bool Sector::Definition::operator==(const Definition& sector_def) const
{
    if (do_stochastic_loss_weighting != sector_def.do_stochastic_loss_weighting)
        return false;
    else if (stochastic_loss_weighting != sector_def.stochastic_loss_weighting)
        return false;
    else if (stopping_decay != sector_def.stopping_decay)
        return false;
    else if (do_continuous_randomization
        != sector_def.do_continuous_randomization)
        return false;
    else if (do_continuous_energy_loss_output
        != sector_def.do_continuous_energy_loss_output)
        return false;
    else if (do_exact_time_calculation != sector_def.do_exact_time_calculation)
        return false;
    else if (only_loss_inside_detector != sector_def.only_loss_inside_detector)
        return false;
    else if (scattering_model != sector_def.scattering_model)
        return false;
    else if (location != sector_def.location)
        return false;
    else if (utility_def != sector_def.utility_def)
        return false;
    else if (cut_settings != sector_def.cut_settings)
        return false;
    else if (*medium_ != *sector_def.medium_)
        return false;
    else if (*geometry_ != *sector_def.geometry_)
        return false;
    return true;
}

bool Sector::Definition::operator!=(const Definition& sector_def) const
{
    return !(*this == sector_def);
}


Sector::Definition::~Definition()
{
}

void Sector::Definition::SetMedium(std::shared_ptr<const Medium> medium)
{
    medium_ = medium;
}

void Sector::Definition::SetGeometry(std::shared_ptr<const Geometry> geometry)
{
    geometry_ = geometry;
}

// ------------------------------------------------------------------------- //
// OStream
// ------------------------------------------------------------------------- //

namespace PROPOSAL{

std::ostream& operator<<(std::ostream& os, Sector const& sector)
{
    std::stringstream ss;
    ss << " Sector (" << &sector << ") ";
    os << Helper::Centered(60, ss.str()) << '\n';

    os << "Sector Definition:\n" << sector.sector_def_ << std::endl;
    os << "Particle Definition:\n" << sector.particle_def_ << std::endl;
    os << "Propagation Utility:\n" << sector.utility_ << std::endl;
    os << "Scattering:\n" << *sector.scattering_ << std::endl;

    os << Helper::Centered(60, "");
    return os;
}

}

// ------------------------------------------------------------------------- //
// Constructors
// ------------------------------------------------------------------------- //

Sector::Sector(const ParticleDef& particle_def, const Definition& sector_def)
    : sector_def_(sector_def)
    , particle_def_(particle_def)
    , utility_(particle_def, sector_def.GetMedium(), sector_def.cut_settings,
          sector_def.utility_def)
    , displacement_calculator_(new UtilityIntegralDisplacement(utility_))
    , interaction_calculator_(new UtilityIntegralInteraction(utility_))
    , decay_calculator_(new UtilityIntegralDecay(utility_))
    , exact_time_calculator_(NULL)
    , cont_rand_(NULL)
    , scattering_(ScatteringFactory::Get().CreateScattering(
          sector_def_.scattering_model, particle_def, utility_))
{
    // These are optional, therfore check NULL
    if (sector_def_.do_exact_time_calculation) {
        exact_time_calculator_ = std::make_shared<UtilityIntegralTime>(utility_);
    }

    if (sector_def_.do_continuous_randomization) {
        cont_rand_ = std::make_shared<ContinuousRandomizer>(utility_);
    }
}

Sector::Sector(const ParticleDef& particle_def, const Definition& sector_def,
    const InterpolationDef& interpolation_def)
    : sector_def_(sector_def)
    , particle_def_(particle_def)
    , utility_(particle_def, sector_def.GetMedium(), sector_def.cut_settings,
          sector_def.utility_def, interpolation_def)
    , displacement_calculator_(
          new UtilityInterpolantDisplacement(utility_, interpolation_def))
    , interaction_calculator_(
          new UtilityInterpolantInteraction(utility_, interpolation_def))
    , decay_calculator_(
          new UtilityInterpolantDecay(utility_, interpolation_def))
    , exact_time_calculator_(NULL)
    , cont_rand_(NULL)
    , scattering_(ScatteringFactory::Get().CreateScattering(
          sector_def_.scattering_model, particle_def, utility_,
          interpolation_def))
{
    // These are optional, therfore check NULL
    if (sector_def_.do_exact_time_calculation) {
        exact_time_calculator_ = std::make_shared<UtilityInterpolantTime>(utility_, interpolation_def);
    }

    if (sector_def_.do_continuous_randomization) {
        cont_rand_ = std::make_shared<ContinuousRandomizer>(utility_, interpolation_def);
    }
}

Sector::Sector(const Sector& sector)
    : sector_def_(sector.sector_def_)
    , particle_def_(sector.particle_def_)
    , utility_(sector.utility_)
    , displacement_calculator_(sector.displacement_calculator_->clone(utility_))
    , interaction_calculator_(sector.interaction_calculator_->clone(utility_))
    , decay_calculator_(sector.decay_calculator_->clone(utility_))
    , exact_time_calculator_(NULL)
    , cont_rand_(sector.cont_rand_)
    , scattering_(sector.scattering_)
{
    // These are optional, therfore check NULL
    if (sector.exact_time_calculator_ != NULL) {
        exact_time_calculator_ = sector.exact_time_calculator_;
    }

}

bool Sector::operator==(const Sector& sector) const
{
    if (sector_def_ != sector.sector_def_)
        return false;
    else if (particle_def_ != sector.particle_def_)
        return false;
    else if (utility_ != sector.utility_)
        return false;
    else if (*cont_rand_ != *sector.cont_rand_)
        return false;
    else if (*scattering_ != *sector.scattering_)
        return false;
    return true;
}

bool Sector::operator!=(const Sector& sector) const
{
    return !(*this == sector);
}

Sector::~Sector()
{

}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %                          Sector Utilities                               %
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double Sector::CalculateTime(const DynamicData& p_condition,
    const double final_energy, const double displacement)
{
    if (exact_time_calculator_) {
        // DensityDistribution Approximation: Use the DensityDistribution at the
        // position of initial energy
        return p_condition.GetTime()
            + exact_time_calculator_->Calculate(
                  p_condition.GetEnergy(), final_energy, 0.0)
            / utility_.GetMedium()->GetDensityDistribution().Evaluate(
                  p_condition.GetPosition());
    }

    return p_condition.GetTime() + displacement / SPEED;
}

void Sector::Scatter(const double displacement, const double initial_energy,
    const double final_energy, Vector3D& position, Vector3D& direction)
{
    if (sector_def_.scattering_model != ScatteringFactory::Enum::NoScattering) {
        Directions directions = scattering_->Scatter(
            displacement, initial_energy, final_energy, position, direction);
        position = position + displacement * directions.u_;
        direction = directions.n_i_;
    } else {
        position = position + displacement * direction;
    }
}

double Sector::ContinuousRandomize(
    const double initial_energy, const double final_energy)
{
    if (cont_rand_) {
        if (final_energy != particle_def_.low) {
            double rnd = RandomGenerator::Get().RandomDouble();
            return cont_rand_->Randomize(initial_energy, final_energy, rnd);
        }
    }
    return final_energy;
}

std::pair<double, int> Sector::MakeStochasticLoss(double particle_energy)
{
    double rnd1 = RandomGenerator::Get().RandomDouble();
    double rnd2 = RandomGenerator::Get().RandomDouble();
    double rnd3 = RandomGenerator::Get().RandomDouble();

    std::pair<double, int> energy_loss;

    // do an optional bias to the randomly sampled energy loss
    if (sector_def_.do_stochastic_loss_weighting) {
        double aux = rnd2
            * std::pow(rnd2, std::abs(sector_def_.stochastic_loss_weighting));
        if (sector_def_.stochastic_loss_weighting > 0) {
            rnd2 = 1 - aux;
        } else {
            rnd2 = aux;
        }
    }

    energy_loss = utility_.StochasticLoss(particle_energy, rnd1, rnd2, rnd3);

    return energy_loss;
}

double Sector::Displacement(const DynamicData& p_condition,
    const double final_energy, const double border_length)
{
    try{
        return displacement_calculator_->Calculate(p_condition.GetEnergy(),
        final_energy, border_length, p_condition.GetPosition(),
        p_condition.GetDirection());
    }
    catch(DensityException& e)
    {
        throw;
    }
    
}

double Sector::BorderLength(const Vector3D& position, const Vector3D& direction)
{
    // loop ueber alle sektoren hoeherer ordnung die im aktuellen sektor liegen
    // und getroffen werden koennten
    return sector_def_.GetGeometry()->DistanceToBorder(position, direction).first;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %                          LossEnergies                                   %
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double Sector::EnergyDecay(const double initial_energy, const double rnd)
{
    double rndd = -std::log(rnd);
    double rnddMin = 0;

    // solving the tracking integral
    if (particle_def_.lifetime < 0) {
        return particle_def_.low;
    }

    rnddMin
        = decay_calculator_->Calculate(initial_energy, particle_def_.low, rndd);

    // evaluating the energy loss
    if (rndd >= rnddMin || rnddMin <= 0) {
        return particle_def_.low;
    }

    return decay_calculator_->GetUpperLimit(initial_energy, rndd);
}

double Sector::EnergyInteraction(const double initial_energy, const double rnd)
{
    double rndi = -std::log(rnd);
    double rndiMin = 0;

    // solving the tracking integral
    rndiMin = interaction_calculator_->Calculate(
        initial_energy, particle_def_.low, rndi);

    if (rndi >= rndiMin || rndiMin <= 0) {
        return particle_def_.low;
    }

    return interaction_calculator_->GetUpperLimit(initial_energy, rndi);
}

double Sector::EnergyMinimal(const double current_energy, const double cut)
{
    return std::min({ current_energy, cut });
}

double Sector::EnergyDistance(
    const double initial_energy, const double distance)
{
    return displacement_calculator_->GetUpperLimit(initial_energy, distance);
}

int Sector::maximizeEnergy(const std::array<double, 4>& LossEnergies)
{
    const auto minmax
        = std::minmax_element(begin(LossEnergies), end(LossEnergies));

    if (*minmax.second == LossEnergies[LossType::Decay]) {
        return LossType::Decay;
    } else {
        return std::distance(LossEnergies.begin(), minmax.second);
    }
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %                               Do Loss                                   %
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

std::shared_ptr<DynamicData> Sector::DoInteraction(
    const DynamicData& p_condition)
{
    std::pair<double, int> stochastic_loss
        = MakeStochasticLoss(p_condition.GetEnergy());

    CrossSection* cross_section
        = utility_.GetCrosssection(stochastic_loss.second);
    std::pair<double, double> deflection_angles
        = cross_section->StochasticDeflection(
            p_condition.GetEnergy(), stochastic_loss.first);

    Vector3D new_direction(p_condition.GetDirection());
    new_direction.deflect(deflection_angles.first, deflection_angles.second);

    double new_energy = p_condition.GetEnergy() - stochastic_loss.first;
    return std::make_shared<DynamicData>(
        stochastic_loss.second,
        p_condition.GetPosition(),
        new_direction,
        new_energy,
        p_condition.GetEnergy(),
        p_condition.GetTime(),
        p_condition.GetPropagatedDistance());
}

std::shared_ptr<DynamicData> Sector::DoDecay(const DynamicData& p_condition)
{

    return std::make_shared<DynamicData>(
        static_cast<int>(InteractionType::Decay), p_condition.GetPosition(),
        p_condition.GetDirection(), p_condition.GetEnergy(),
        p_condition.GetParentParticleEnergy(), p_condition.GetTime(),
        p_condition.GetPropagatedDistance());
}

std::shared_ptr<DynamicData> Sector::DoContinuous(
    const DynamicData& p_condition, double final_energy, double displacement)
{

    double initial_energy{ p_condition.GetEnergy() };

    double dist = p_condition.GetPropagatedDistance() + displacement;
    double time = CalculateTime(p_condition, final_energy, displacement);

    Vector3D position{ p_condition.GetPosition() };
    Vector3D direction{ p_condition.GetDirection() };

    Scatter(displacement, p_condition.GetEnergy(), final_energy, position,
        direction);
    final_energy = ContinuousRandomize(p_condition.GetEnergy(), final_energy);

    return std::make_shared<DynamicData>(
        static_cast<int>(InteractionType::ContinuousEnergyLoss), position,
        direction, final_energy, initial_energy, time, dist);
}

Secondaries Sector::Propagate(
    const DynamicData& p_initial, double border_distance, const double minimal_energy)
{
    Secondaries secondaries(std::make_shared<ParticleDef>(particle_def_));

    auto p_condition = std::make_shared<DynamicData>(p_initial);
    double dist_limit{ p_initial.GetPropagatedDistance() + border_distance };
    double rnd;
    int minimalLoss;
    std::array<double, 4> LossEnergies;
    double displacement;

    while (true) {
        rnd = RandomGenerator::Get().RandomDouble();
        LossEnergies[LossType::Decay]
            = EnergyDecay(p_condition->GetEnergy(), rnd);

        rnd = RandomGenerator::Get().RandomDouble();
        LossEnergies[LossType::Interaction]
            = EnergyInteraction(p_condition->GetEnergy(), rnd);

        border_distance = dist_limit - p_condition->GetPropagatedDistance();
        LossEnergies[LossType::Distance]
            = EnergyDistance(p_condition->GetEnergy(), border_distance);

        LossEnergies[LossType::MinimalE]
            = EnergyMinimal(p_condition->GetEnergy(), minimal_energy);

        minimalLoss = maximizeEnergy(LossEnergies);

        if (minimalLoss == LossType::Distance)
        {
            displacement = border_distance;
        }
        else
        {
            try{
                displacement = Displacement(*p_condition, LossEnergies[minimalLoss], border_distance);
            }
            catch(DensityException& e){
                // due to numerical instabilities
                // the calculated Energy to the distance to the sector border
                // may be slightly bigger than an interaction length, which could
                // exceed the sector border.
                minimalLoss = LossType::Distance;
                displacement = border_distance;
            }

            if(std::abs(displacement - border_distance) < PARTICLE_POSITION_RESOLUTION){
                minimalLoss = LossType::Distance;
                displacement = border_distance;
            }
        }

        p_condition
            = DoContinuous(*p_condition, LossEnergies[minimalLoss], displacement);
        if (sector_def_.do_continuous_energy_loss_output)
            secondaries.push_back(*p_condition);

        if (minimalLoss == LossType::Interaction)
        {
            p_condition = DoInteraction(*p_condition);
            secondaries.push_back(*p_condition);
        }
        else
        {
            break;
        }
    };

    if (minimalLoss == LossType::Decay)
    {
        p_condition = DoDecay(*p_condition);
    }

    secondaries.push_back(*p_condition);

    return secondaries;
}
