#include "PROPOSAL/crosssection/Factories/PhotoproductionFactory.h"
#include "PROPOSAL/crosssection/parametrization/Photoproduction.h"
#include "PROPOSAL/crosssection/CrossSectionBuilder.h"
#include "PROPOSAL/crosssection/CrossSectionMultiplier.h"
#include <nlohmann/json.hpp>

using namespace PROPOSAL;
using photo_func_ptr = cross_ptr (*)(const ParticleDef&, const Medium&);

template <typename Param>
cross_ptr create_photoproduction(const ParticleDef& p_def, const Medium& medium)
{
    auto param = Param();
    return make_crosssection(param, p_def, medium, nullptr, false);
}

template <typename Param>
auto init_param = std::make_pair(std::string(crosssection::ParametrizationName<Param>::value), create_photoproduction<Param>);

std::map<std::string, photo_func_ptr, Helper::case_insensitive_comp> photoproduction_map = {
        init_param<crosssection::PhotoproductionZeus>,
        init_param<crosssection::PhotoproductionBezrukovBugaev>,
        init_param<crosssection::PhotoproductionCaldwell>,
        init_param<crosssection::PhotoproductionKokoulin>,
        init_param<crosssection::PhotoproductionRhode>
};

namespace PROPOSAL {
    cross_ptr make_photoproduction(const ParticleDef &p_def, const Medium &medium,
                                   const nlohmann::json &config) {
        if (!config.contains("parametrization"))
            throw std::logic_error("No parametrization passed for photoproduction");
        std::string param_name = config["parametrization"];

        auto it = photoproduction_map.find(param_name);
        if (it == photoproduction_map.end())
            throw std::logic_error("Unknown parametrization for photoproduction");

        auto cross = it->second(p_def, medium);

        double multiplier = config.value("multiplier", 1.0);
        if (multiplier != 1.0)
            return make_crosssection_multiplier(std::shared_ptr<CrossSectionBase>(std::move(cross)), multiplier);
        return cross;
    }
}
