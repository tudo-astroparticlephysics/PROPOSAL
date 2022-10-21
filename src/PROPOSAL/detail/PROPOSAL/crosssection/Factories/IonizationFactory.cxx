#include "PROPOSAL/crosssection/Factories/IonizationFactory.h"
#include "PROPOSAL/crosssection/parametrization/Ionization.h"
#include "PROPOSAL/crosssection/CrossSectionBuilder.h"
#include "PROPOSAL/crosssection/CrossSectionMultiplier.h"
#include <nlohmann/json.hpp>

using namespace PROPOSAL;
using ioniz_func_ptr = cross_ptr (*)(const ParticleDef&, const Medium&, std::shared_ptr<const EnergyCutSettings>, bool, bool);

template <typename Param>
cross_ptr create_ionization(const ParticleDef& p_def, const Medium& medium,
                            std::shared_ptr<const EnergyCutSettings> cuts,
                            bool interpol, bool enable_density_correction)
{
    auto param = Param(*cuts, enable_density_correction);
    return make_crosssection(param, p_def, medium, cuts, interpol);
}

template <typename Param>
auto init_param = std::make_pair(std::string(crosssection::ParametrizationName<Param>::value), create_ionization<Param>);

std::map<std::string, ioniz_func_ptr, Helper::case_insensitive_comp> ioniz_map = {
        init_param<crosssection::IonizBetheBlochRossi>,
        init_param<crosssection::IonizBergerSeltzerBhabha>,
        init_param<crosssection::IonizBergerSeltzerMoller>
};

namespace PROPOSAL {
    cross_ptr make_ionization(const ParticleDef &p_def, const Medium &medium,
                             std::shared_ptr<const EnergyCutSettings> cuts,
                             bool interpol, const nlohmann::json &config) {
        if (!config.contains("parametrization"))
            throw std::logic_error("No parametrization passed for ionization");
        std::string param_name = config["parametrization"];

        auto it = ioniz_map.find(param_name);
        if (it == ioniz_map.end())
            throw std::logic_error("Unknown parametrization for ionization");
        bool enable_density_correction = config.value("enable_density_correction", true);
        auto cross = it->second(p_def, medium, cuts, interpol, enable_density_correction);

        double multiplier = config.value("multiplier", 1.0);
        if (multiplier != 1.0)
            return make_crosssection_multiplier(std::shared_ptr<CrossSectionBase>(std::move(cross)), multiplier);
        return cross;
    }
}
