#include "PROPOSAL/crosssection/Factories/ComptonFactory.h"
#include "PROPOSAL/crosssection/parametrization/Compton.h"
#include "PROPOSAL/crosssection/CrossSectionBuilder.h"
#include "PROPOSAL/crosssection/CrossSectionMultiplier.h"

using namespace PROPOSAL;
using compton_func_ptr = cross_ptr (*)(const ParticleDef&, const Medium&, std::shared_ptr<const EnergyCutSettings>, bool);

template <typename Param>
cross_ptr create_compton(const ParticleDef& p_def, const Medium& medium,
                         std::shared_ptr<const EnergyCutSettings> cuts,
                         bool interpol)
{
    auto param = Param();
    return make_crosssection(param, p_def, medium, cuts, interpol);
}

template <typename Param>
auto init_param = std::make_pair(std::string(crosssection::ParametrizationName<Param>::value), create_compton<Param>);

std::map<std::string, compton_func_ptr, Helper::case_insensitive_comp> compton_map = {
        init_param<crosssection::ComptonKleinNishina>
};

namespace PROPOSAL {
    cross_ptr make_compton(const ParticleDef &p_def, const Medium &medium,
                           std::shared_ptr<const EnergyCutSettings> cuts,
                           bool interpol, const nlohmann::json &config) {
        if (!config.contains("parametrization"))
            throw std::logic_error("No parametrization passed for compton");
        std::string param_name = config["parametrization"];

        auto it = compton_map.find(param_name);
        if (it == compton_map.end())
            throw std::logic_error("Unknown parametrization for compton");

        auto cross = it->second(p_def, medium, cuts, interpol);

        double multiplier = config.value("multiplier", 1.0);
        if (multiplier != 1.0)
            return make_crosssection_multiplier(std::shared_ptr<CrossSectionBase>(std::move(cross)), multiplier);
        return cross;
    }
}
