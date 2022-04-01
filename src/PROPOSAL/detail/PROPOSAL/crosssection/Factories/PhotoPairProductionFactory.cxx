#include "PROPOSAL/crosssection/Factories/PhotoPairProductionFactory.h"
#include "PROPOSAL/crosssection/parametrization/PhotoPairProduction.h"
#include "PROPOSAL/crosssection/CrossSectionBuilder.h"
#include "PROPOSAL/crosssection/CrossSectionMultiplier.h"
#include <nlohmann/json.hpp>

using namespace PROPOSAL;
using photopair_func_ptr = cross_ptr (*)(const ParticleDef&, const Medium&,
        bool, bool, double);

template <typename Param>
cross_ptr create_photopairproduction(
        const ParticleDef& p_def, const Medium& medium, bool interpol,
        bool lpm, double density_correction)
{
    auto param = Param(lpm, p_def, medium, density_correction);
    return make_crosssection(param, p_def, medium, nullptr, interpol);
}

template <typename Param>
auto init_param = std::make_pair(std::string(crosssection::ParametrizationName<Param>::value), create_photopairproduction<Param>);

std::map<std::string, photopair_func_ptr, Helper::case_insensitive_comp> photopair_map = {
        init_param<crosssection::PhotoPairTsai>,
        init_param<crosssection::PhotoPairKochMotz>
};

namespace PROPOSAL {
    cross_ptr make_photopairproduction(const ParticleDef &p_def, const Medium &medium,
                                       bool interpol, const nlohmann::json &config,
                                       double density_correction) {
        if (!config.contains("parametrization"))
            throw std::logic_error("No parametrization passed for photopairproduction");
        std::string param_name = config["parametrization"];
        bool lpm = config.value("lpm", true);

        auto it = photopair_map.find(param_name);
        if (it == photopair_map.end())
            throw std::logic_error("Unknown parametrization for photopairproduction");

        auto cross = it->second(p_def, medium, interpol, lpm, density_correction);

        double multiplier = config.value("multiplier", 1.0);
        if (multiplier != 1.0)
            return make_crosssection_multiplier(std::shared_ptr<CrossSectionBase>(std::move(cross)), multiplier);
        return cross;
    }
}
