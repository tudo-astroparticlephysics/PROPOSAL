#include "PROPOSAL/crosssection/Factories/PhotoeffectFactory.h"
#include "PROPOSAL/crosssection/parametrization/Photoeffect.h"
#include "PROPOSAL/crosssection/CrossSectionBuilder.h"
#include "PROPOSAL/crosssection/CrossSectionMultiplier.h"
#include <nlohmann/json.hpp>

using namespace PROPOSAL;
using photoeffect_func_ptr = cross_ptr (*)(const ParticleDef&, const Medium&);

template <typename Param>
cross_ptr create_photoeffect(const ParticleDef& p_def, const Medium& medium)
{
    auto param = Param();
    return make_crosssection(param, p_def, medium, nullptr, false);
}

template <typename Param>
auto init_param = std::make_pair(std::string(crosssection::ParametrizationName<Param>::value), create_photoeffect<Param>);

std::map<std::string, photoeffect_func_ptr, Helper::case_insensitive_comp> photoeffect_map = {
        init_param<crosssection::PhotoeffectSauter>
};

namespace PROPOSAL {
    cross_ptr make_photoeffect(const ParticleDef &p_def, const Medium &medium,
                               const nlohmann::json &config) {
        if (!config.contains("parametrization"))
            throw std::logic_error("No parametrization passed for photoeffect");
        std::string param_name = config["parametrization"];

        auto it = photoeffect_map.find(param_name);
        if (it == photoeffect_map.end())
            throw std::logic_error("Unknown parametrization for photoeffect");

        auto cross = it->second(p_def, medium);

        double multiplier = config.value("multiplier", 1.0);
        if (multiplier != 1.0)
            return make_crosssection_multiplier(std::shared_ptr<CrossSectionBase>(std::move(cross)), multiplier);
        return cross;
    }
}
