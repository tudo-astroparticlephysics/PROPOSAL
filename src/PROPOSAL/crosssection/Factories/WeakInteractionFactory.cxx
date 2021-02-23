#include "PROPOSAL/crosssection/Factories/WeakInteractionFactory.h"
#include "PROPOSAL/crosssection/parametrization/WeakInteraction.h"
#include "PROPOSAL/crosssection/CrossSectionBuilder.h"
#include "PROPOSAL/crosssection/CrossSectionMultiplier.h"

using namespace PROPOSAL;
using weak_func_ptr = cross_ptr (*)(const ParticleDef&, const Medium&, bool);

template <typename Param>
cross_ptr create_weakinteraction(const ParticleDef& p_def, const Medium& medium, bool interpol)
{
    auto param = Param();
    return make_crosssection(param, p_def, medium, nullptr, interpol);
}

template <typename Param>
auto init_param = std::make_pair(std::string(crosssection::ParametrizationName<Param>::value), create_weakinteraction<Param>);

std::map<std::string, weak_func_ptr, Helper::case_insensitive_comp> weak_map = {
        init_param<crosssection::WeakCooperSarkarMertsch>
};

namespace PROPOSAL {
    cross_ptr make_weakinteraction(const ParticleDef &p_def, const Medium &medium,
                                   bool interpol, const nlohmann::json &config) {
        if (!config.contains("parametrization"))
            throw std::logic_error("No parametrization passed for weak interaction");
        std::string param_name = config["parametrization"];

        auto it = weak_map.find(param_name);
        if (it == weak_map.end())
            throw std::logic_error("Unknown parametrization for weak interaction");

        auto cross = it->second(p_def, medium, interpol);

        double multiplier = config.value("multiplier", 1.0);
        if (multiplier != 1.0)
            return make_crosssection_multiplier(std::shared_ptr<CrossSectionBase>(std::move(cross)), multiplier);
        return cross;
    }
}
