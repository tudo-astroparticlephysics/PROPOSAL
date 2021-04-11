#pragma once

#include "PROPOSAL/Secondaries.h"
#include <nlohmann/json.hpp>
#include <unordered_map>

namespace PROPOSAL {
struct ParticleDef;

class Geometry;

using Sector = std::tuple<std::shared_ptr<const Geometry>, PropagationUtility,
    std::shared_ptr<const Density_distr>>;

struct CrossSectionBase;
}

namespace PROPOSAL {
class Propagator {
public:
    Propagator(const ParticleDef& p_def, const nlohmann::json& config);
    Propagator(const ParticleDef& p_def, const std::string& config_file)
        : Propagator(p_def, ParseConfig(config_file))
    {
    }
    Propagator(const ParticleDef&, std::vector<Sector> sectors);

    Secondaries Propagate(const ParticleState& initial_particle,
        double max_distance = 1e20, double min_energy = 0.,
        unsigned int hierarchy_condition = 0);
    enum { GEOMETRY, UTILITY, DENSITY_DISTR };

private:
    InteractionType DoStochasticInteraction(
        ParticleState&, PropagationUtility&, std::function<double()>);
    int AdvanceParticle(ParticleState& p_cond, double E_f, double max_distance,
        std::function<double()> rnd, Sector& current_sector);
    double CalculateDistanceToBorder(const Vector3D& particle_position,
        const Vector3D& particle_direction, const Geometry& current_geometry);
    int maximize(const std::array<double, 3>& InteractionEnergies);
    int minimize(const std::array<double, 3>& AdvanceDistances);
    Sector GetCurrentSector(
        const Vector3D& particle_position, const Vector3D& particle_direction);
    // Global settings
    struct GlobalSettings {
        GlobalSettings();
        GlobalSettings(const nlohmann::json& config_global);
        nlohmann::json cross;
        nlohmann::json scattering;
        std::shared_ptr<EnergyCutSettings> cuts = nullptr;
        std::shared_ptr<Medium> medium = nullptr;
        bool do_exact_time;
        bool do_interpolation;
    };

    // Initializing methods
    static nlohmann::json ParseConfig(const std::string& config_file);
    void InitializeSectorFromJSON(
        const ParticleDef&, const nlohmann::json&, GlobalSettings);

    PropagationUtility::Collection CreateUtility(
        std::vector<std::shared_ptr<CrossSectionBase>> crosss,
        std::shared_ptr<Medium> medium, bool do_cont_rand, bool do_interpol,
        bool do_exact_time, nlohmann::json scatter);

    std::vector<std::shared_ptr<CrossSectionBase>> CreateCrossSectionList(
        const ParticleDef& p_def, const Medium& medium,
        std::shared_ptr<const EnergyCutSettings> cuts, bool interpolate,
        double density_correction, const nlohmann::json& config);

    ParticleDef p_def;
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
