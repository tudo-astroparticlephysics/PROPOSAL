#include "PROPOSAL/crosssection/Factories/EpairProductionFactory.h"
#include "PROPOSAL/crosssection/parametrization/EpairProduction.h"
#include "PROPOSAL/crosssection/CrossSectionBuilder.h"
#include "PROPOSAL/crosssection/CrossSectionMultiplier.h"

using namespace PROPOSAL;
using epair_func_ptr = cross_ptr (*)(const ParticleDef&, const Medium&,
                                     std::shared_ptr<const EnergyCutSettings>, bool, bool, double);

template <typename Param>
cross_ptr create_epair(const ParticleDef& p_def, const Medium& medium,
                       std::shared_ptr<const EnergyCutSettings> cuts, bool lpm,
                       bool interpol, double density_correction)
{
    auto param = Param(lpm, p_def, medium, density_correction);
    return make_crosssection(param, p_def, medium, cuts, interpol);
}

template <typename Param>
auto init_param = std::make_pair(std::string(crosssection::ParametrizationName<Param>::value), create_epair<Param>);

std::map<std::string, epair_func_ptr, Helper::case_insensitive_comp> epair_map = {
        init_param<crosssection::EpairKelnerKokoulinPetrukhin>,
        init_param<crosssection::EpairSandrockSoedingreksoRhode>,
        init_param<crosssection::EpairForElectronPositron>
};

namespace PROPOSAL {
    cross_ptr make_epairproduction(const ParticleDef &p_def, const Medium &medium,
                                  std::shared_ptr<const EnergyCutSettings> cuts,
                                  bool interpol, const nlohmann::json& config,
                                  double density_correction) {
        if (!config.contains("parametrization"))
            throw std::logic_error("No parametrization passed for epairproduction");

        std::string param_name = config["parametrization"];
        bool lpm = config.value("lpm", true);

        auto it = epair_map.find(param_name);
        if (it == epair_map.end())
            throw std::logic_error("Unknown parametrization for epairproduction");

        auto cross = it->second(p_def, medium, cuts, lpm, interpol, density_correction);

        double multiplier = config.value("multiplier", 1.0);
        if (multiplier != 1.0)
            return make_crosssection_multiplier(std::shared_ptr<CrossSectionBase>(std::move(cross)), multiplier);
        return cross;
    }
}
