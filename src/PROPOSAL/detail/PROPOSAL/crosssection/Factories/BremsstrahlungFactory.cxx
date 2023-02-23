#include "PROPOSAL/crosssection/Factories/BremsstrahlungFactory.h"
#include "PROPOSAL/crosssection/parametrization/Bremsstrahlung.h"
#include "PROPOSAL/crosssection/CrossSectionBuilder.h"
#include "PROPOSAL/crosssection/CrossSectionMultiplier.h"
#include <nlohmann/json.hpp>

using namespace PROPOSAL;
using brems_func_ptr = cross_ptr (*)(const ParticleDef&, const Medium&,
        std::shared_ptr<const EnergyCutSettings>, bool, bool, double, bool);

template <typename Param>
cross_ptr create_brems(const ParticleDef& p_def, const Medium& medium,
                       std::shared_ptr<const EnergyCutSettings> cuts, bool lpm,
                       bool interpol, double density_correction, bool TM_effect)
{
    auto param = Param(lpm, p_def, medium, density_correction, TM_effect);
    return make_crosssection(param, p_def, medium, cuts, interpol);
}

template <typename Param>
auto init_param = std::make_pair(std::string(crosssection::ParametrizationName<Param>::value), create_brems<Param>);

std::map<std::string, brems_func_ptr, Helper::case_insensitive_comp> brems_map = {
        init_param<crosssection::BremsKelnerKokoulinPetrukhin>,
        init_param<crosssection::BremsPetrukhinShestakov>,
        init_param<crosssection::BremsCompleteScreening>,
        init_param<crosssection::BremsAndreevBezrukovBugaev>,
        init_param<crosssection::BremsSandrockSoedingreksoRhode>,
        init_param<crosssection::BremsElectronScreening>,
};

namespace PROPOSAL {
    cross_ptr make_bremsstrahlung(const ParticleDef &p_def, const Medium &medium,
                                  std::shared_ptr<const EnergyCutSettings> cuts,
                                  bool interpol, const nlohmann::json& config,
                                  double density_correction) {
        if (!config.contains("parametrization"))
            throw std::logic_error("No parametrization passed for bremsstrahlung");

        std::string param_name = config["parametrization"];
        bool lpm = config.value("lpm", true);
        bool TM_effect = config.value("tm_effect", true);
        auto it = brems_map.find(param_name);
        if (it == brems_map.end())
            throw std::logic_error("Unknown parametrization for bremsstrahlung");

        auto cross = it->second(p_def, medium, cuts, lpm, interpol, density_correction, TM_effect);

        double multiplier = config.value("multiplier", 1.0);
        if (multiplier != 1.0)
            return make_crosssection_multiplier(std::shared_ptr<CrossSectionBase>(std::move(cross)), multiplier);
        return cross;
        }
    }
