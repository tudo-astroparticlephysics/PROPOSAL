#pragma once

#include <unordered_map>

#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/density_distr/density_distr.h"
#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/json.hpp"
#include "PROPOSAL/geometry/Geometry.h"
#include "PROPOSAL/propagation_utility/TimeBuilder.h"
#include "PROPOSAL/geometry/GeometryFactory.h"
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
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/scattering/multiple_scattering/ScatteringFactory.h"
#include "PROPOSAL/Secondaries.h"

using std::get;
using std::make_shared;
namespace PROPOSAL {

using Sector = std::tuple<std::shared_ptr<const Geometry>, PropagationUtility,
    std::shared_ptr<const Density_distr>>;

/* class Secondaries; */

class Propagator
{
public:
    template <typename P>
    Propagator(P&& p_def, const nlohmann::json& config) : p_def(p_def)
    {
        GlobalSettings global;
        if (config.contains("global"))
            global = GlobalSettings(config["global"]);
        if (config.contains("interpolation")) {
            interpol_def_global
                    = std::make_shared<InterpolationDef>(config["interpolation"]);
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
    template <typename P>
    Propagator(P&& p_def, const std::string& config_file)
    : Propagator(p_def, ParseConfig(config_file)) {}
    Propagator(ParticleDef const&, std::vector<Sector> sectors);

    Secondaries Propagate(const ParticleState& initial_particle,
            double max_distance = 1e20, double min_energy = 0.);
    enum {GEOMETRY, UTILITY, DENSITY_DISTR};

private:
    InteractionType DoStochasticInteraction(ParticleState&, PropagationUtility&,
                                            std::function<double()>);
    int AdvanceParticle(ParticleState& p_cond, double E_f, double max_distance,
                        std::function<double()> rnd, Sector& current_sector);
    double CalculateDistanceToBorder(const Vector3D& particle_position,
            const Vector3D& particle_direction,
            const Geometry& current_geometry);
    int maximize(const std::array<double, 3>& InteractionEnergies);
    int minimize(const std::array<double, 3>& AdvanceDistances);
    Sector GetCurrentSector(const Vector3D& particle_position,
                            const Vector3D& particle_direction);
    //Global settings
    struct GlobalSettings{
        GlobalSettings();
        GlobalSettings(const nlohmann::json& config_global);
        nlohmann::json cross;
        std::string scattering;
        std::shared_ptr<EnergyCutSettings> cuts = nullptr;
        std::shared_ptr<Medium> medium = nullptr;
        bool do_exact_time;
        bool do_interpolation;
    };

    //Initializing methods
    static nlohmann::json ParseConfig(const std::string& config_file);
    template <typename P>
    void InitializeSectorFromJSON(P&& p_def, const nlohmann::json& json_sector,
                                  GlobalSettings global)
    {
        bool do_interpolation
                = json_sector.value("do_interpolation", global.do_interpolation);
        bool do_exact_time = json_sector.value("exact_time", global.do_exact_time);
        std::string scattering = json_sector.value("scattering", global.scattering);
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
        nlohmann::json density_distr = {{"type", "homogeneous"},
                                        {"mass_density", medium->GetMassDensity()}};
        if(json_sector.contains("density_distribution"))
            density_distr = json_sector["density_distribution"];

        auto cross_config = json_sector.value("CrossSections", global.cross);
        PropagationUtility::Collection collection;
        if (!cross_config.empty()) {
            double density_correction = density_distr.value("mass_density", medium->GetMassDensity());
            density_correction /= medium->GetMassDensity();
            auto crosss = CreateCrossSectionList(p_def, *medium, cuts,
                                                 do_interpolation, density_correction,
                                                 cross_config);
            collection = CreateUtility(crosss, medium, cuts->GetContRand(),
                                       do_interpolation, do_exact_time, scattering);
        } else {
            auto std_crosss = GetStdCrossSections(p_def, *medium, cuts,
                                                  do_interpolation);
            collection = CreateUtility(std_crosss, medium, cuts->GetContRand(),
                                       do_interpolation, do_exact_time, scattering);
        }
        auto utility = PropagationUtility(collection);

        if (json_sector.contains("geometries")) {
            assert(json_sector["geometries"].is_array());
            for (const auto& json_geometry : json_sector.at("geometries")) {
                auto geometry = CreateGeometry(json_geometry);
                auto density = CreateDensityDistribution(density_distr);
                sector_list.emplace_back(std::make_tuple(geometry, utility, density));
            }
        } else {
            throw std::invalid_argument(
                    "At least one geometry must be defined for each sector");
        }
    }

    template<typename CrossVec>
    PropagationUtility::Collection CreateUtility(
            CrossVec&& crosss, std::shared_ptr<Medium> medium, bool do_cont_rand,
            bool do_interpol, bool do_exact_time, std::string scatter)
    {
        PropagationUtility::Collection def;
        def.interaction_calc = make_interaction(crosss, do_interpol);
        def.displacement_calc = make_displacement(crosss, do_interpol);
        def.scattering = make_scattering(scatter, p_def, *medium, crosss,
                                         do_interpol);
        if (std::isfinite(p_def.lifetime))
            def.decay_calc = make_decay(crosss, p_def, do_interpol);
        if (do_cont_rand)
            def.cont_rand = make_contrand(crosss, do_interpol);
        if(do_exact_time) {
            def.time_calc = make_time(crosss, p_def, do_interpol);
        } else {
            def.time_calc = make_shared<ApproximateTimeBuilder>();
        }
        return def;
    }

    template <typename P, typename M>
    crosssection_list_t<P, M> CreateCrossSectionList(
            P&& p_def, M&& medium, shared_ptr<const EnergyCutSettings> cuts,
            bool interpolate, double density_correction, const nlohmann::json& config) {
        crosssection_list_t<P, M> cross;

        if (config.contains("annihilation"))
            cross.emplace_back(crosssection::make_annihilation(p_def, medium,
                                                               interpolate, config["annihilation"]));
        if (config.contains("brems"))
            cross.emplace_back(crosssection::make_bremsstrahlung(p_def, medium,
                                                                 cuts, interpolate, config["brems"],
                                                                 density_correction));
        if (config.contains("compton"))
            cross.emplace_back(crosssection::make_compton(p_def, medium,
                                                          cuts, interpolate, config["compton"]));
        if (config.contains("epair"))
            cross.emplace_back(crosssection::make_epairproduction(p_def, medium,
                                                                  cuts, interpolate, config["epair"],
                                                                  density_correction));
        if (config.contains("ioniz"))
            cross.emplace_back(crosssection::make_ionization(p_def, medium,
                                                             cuts, interpolate, config["ioniz"]));
        if (config.contains("mupair"))
            cross.emplace_back(crosssection::make_mupairproduction(p_def, medium,
                                                                   cuts, interpolate, config["mupair"]));
        if (config.contains("photo")) {
            try { cross.emplace_back(crosssection::make_photonuclearreal(p_def, medium,
                                                                         cuts, interpolate, config["photo"]));
            } catch (std::invalid_argument &e) {
                cross.emplace_back(crosssection::make_photonuclearQ2(p_def, medium,
                                                                     cuts, interpolate, config["photo"]));
            }
        }
        if (config.contains("photopair"))
            cross.emplace_back(crosssection::make_photopairproduction(p_def, medium,
                                                                      interpolate, config["photopair"]));
        if(config.contains("weak"))
            cross.emplace_back(crosssection::make_weakinteraction(p_def, medium,
                                                                  interpolate, config["weak"]));

        return cross;
    }

    ParticleDef p_def;
    std::shared_ptr<InterpolationDef> interpol_def_global = nullptr;
    enum Type : int {
        MinimalE = 0,
        Decay = 1,
        Stochastic = 2,
    };
    enum AdvancementType : int {
        ReachedInteraction = 0,
        ReachedMaxDistance = 1,
        ReachedBorder = 2
    };

    std::vector<Sector> sector_list;

};

} // namespace PROPOSAL
