#include "PROPOSAL/crosssection/Factories/MupairProductionFactory.h"
#include "PROPOSAL/crosssection/parametrization/MupairProduction.h"
#include "PROPOSAL/crosssection/CrossSectionBuilder.h"
#include "PROPOSAL/crosssection/CrossSectionMultiplier.h"

using namespace PROPOSAL;
using mupair_func_ptr = cross_ptr (*)(const ParticleDef&, const Medium&, std::shared_ptr<const EnergyCutSettings>, bool);

template <typename Param>
cross_ptr create_mupairproduction(const ParticleDef& p_def, const Medium& medium,
                                  std::shared_ptr<const EnergyCutSettings> cuts,
                                  bool interpol)
{
    auto param = Param();
    return make_crosssection(param, p_def, medium, cuts, interpol);
}

template <typename Param>
auto init_param = std::make_pair(std::string(crosssection::ParametrizationName<Param>::value), create_mupairproduction<Param>);

std::map<std::string, mupair_func_ptr, Helper::case_insensitive_comp> mupair_map = {
        init_param<crosssection::MupairKelnerKokoulinPetrukhin>
};

namespace PROPOSAL {
    cross_ptr make_mupairproduction(const ParticleDef &p_def, const Medium &medium,
                                    std::shared_ptr<const EnergyCutSettings> cuts,
                                    bool interpol, const nlohmann::json &config) {
        if (!config.contains("parametrization"))
            throw std::logic_error("No parametrization passed for mupairproduction");

        std::string param_name = config["parametrization"];

        auto it = mupair_map.find(param_name);
        if (it == mupair_map.end())
            throw std::logic_error("Unknown parametrization for mupairproduction");

        auto cross = it->second(p_def, medium, cuts, interpol);

        double multiplier = config.value("multiplier", 1.0);
        if (multiplier != 1.0)
            return make_crosssection_multiplier(std::shared_ptr<CrossSectionBase>(std::move(cross)), multiplier);
        return cross;
    }
}
