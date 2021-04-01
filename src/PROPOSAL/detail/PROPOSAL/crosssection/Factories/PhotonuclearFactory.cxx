#include "PROPOSAL/crosssection/Factories/PhotonuclearFactory.h"
#include "PROPOSAL/crosssection/parametrization/PhotoQ2Integration.h"
#include "PROPOSAL/crosssection/parametrization/PhotoRealPhotonAssumption.h"
#include "PROPOSAL/crosssection/CrossSectionBuilder.h"
#include "PROPOSAL/crosssection/CrossSectionMultiplier.h"

using namespace PROPOSAL;
using photoQ2_func_ptr = cross_ptr (*)(const ParticleDef&, const Medium&, std::shared_ptr<const EnergyCutSettings>, std::shared_ptr<crosssection::ShadowEffect>, bool);
using photoreal_func_ptr = cross_ptr (*)(const ParticleDef&, const Medium&, std::shared_ptr<const EnergyCutSettings>, bool, bool);
using shadow_func_ptr = std::shared_ptr<crosssection::ShadowEffect> (*)();

template <typename Param>
cross_ptr create_photonuclearQ2(const ParticleDef& p_def, const Medium& medium,
                                std::shared_ptr<const EnergyCutSettings> cuts,
                                std::shared_ptr<crosssection::ShadowEffect> shadow,
                                bool interpol)
{
    auto param = Param(shadow);
    return make_crosssection(param, p_def, medium, cuts, interpol);
}

template <typename Param>
cross_ptr create_photoreal(const ParticleDef& p_def, const Medium& medium,
                           std::shared_ptr<const EnergyCutSettings> cuts,
                           bool hard_component, bool interpol)
{
    auto param = Param(hard_component);
    return make_crosssection(param, p_def, medium, cuts, interpol);
}

template <typename Param>
auto init_param_Q2 = std::make_pair(std::string(crosssection::ParametrizationName<Param>::value), create_photonuclearQ2<Param>);

std::map<std::string, photoQ2_func_ptr, Helper::case_insensitive_comp> photoQ2_map = {
    init_param_Q2<crosssection::PhotoAbramowiczLevinLevyMaor91>,
    init_param_Q2<crosssection::PhotoAbramowiczLevinLevyMaor97>,
    init_param_Q2<crosssection::PhotoButkevichMikheyev>,
    init_param_Q2<crosssection::PhotoRenoSarcevicSu>,
    init_param_Q2<crosssection::PhotoAbtFT>,
    init_param_Q2<crosssection::PhotoBlockDurandHa>
};

template <typename Param>
auto init_param_photoreal = std::make_pair(std::string(crosssection::ParametrizationName<Param>::value), create_photoreal<Param>);

std::map<std::string, photoreal_func_ptr, Helper::case_insensitive_comp> photoreal_map = {
        init_param_photoreal<crosssection::PhotoZeus>,
        init_param_photoreal<crosssection::PhotoBezrukovBugaev>,
        init_param_photoreal<crosssection::PhotoKokoulin>,
        init_param_photoreal<crosssection::PhotoRhode>,
};

template <typename Shadow>
std::shared_ptr<crosssection::ShadowEffect> create_shadow() {
    return std::make_shared<Shadow>();
}

static std::map<std::string, shadow_func_ptr, Helper::case_insensitive_comp> shadow_map = {
    { "duttarenosarcevicseckel", create_shadow<crosssection::ShadowDuttaRenoSarcevicSeckel>},
    { "butkevichmikheyev", create_shadow<crosssection::ShadowButkevichMikheyev>}
};

namespace PROPOSAL {
    cross_ptr make_photonuclearQ2(const ParticleDef &p_def, const Medium &medium,
                                  std::shared_ptr<const EnergyCutSettings> cuts,
                                  bool interpol, const nlohmann::json &config) {
        if (!config.contains("parametrization"))
            throw std::logic_error("No parametrization passed for photonuclear");

        std::string param_name = config["parametrization"];
        std::string shadow_name = config.value("shadow", "ButkevichMikheyev");


        auto it = photoQ2_map.find(param_name);
        if (it == photoQ2_map.end())
            throw std::invalid_argument("Unknown parametrization for photonuclear");

        auto it_shadow = shadow_map.find(shadow_name);
        if (it_shadow == shadow_map.end())
            throw std::logic_error("Shadow effect name unknown");
        auto cross = it->second(p_def, medium, cuts, it_shadow->second(), interpol);

        double multiplier = config.value("multiplier", 1.0);
        if (multiplier != 1.0)
            return make_crosssection_multiplier(std::shared_ptr<CrossSectionBase>(std::move(cross)), multiplier);
        return cross;
    }
}

namespace PROPOSAL {
    cross_ptr make_photonuclearreal(const ParticleDef &p_def, const Medium &medium,
                                    std::shared_ptr<const EnergyCutSettings> cuts,
                                    bool interpol, const nlohmann::json &config) {
        if (!config.contains("parametrization"))
            throw std::logic_error("No parametrization passed for photonuclear");

        std::string param_name = config["parametrization"];
        bool hard_component = config.value("hard_component", true);

        auto it = photoreal_map.find(param_name);
        if (it == photoreal_map.end())
            throw std::invalid_argument("Unknown parametrization for photonuclear");

        auto cross = it->second(p_def, medium, cuts, hard_component, interpol);

        double multiplier = config.value("multiplier", 1.0);
        if (multiplier != 1.0)
            return make_crosssection_multiplier(std::shared_ptr<CrossSectionBase>(std::move(cross)), multiplier);
        return cross;
    }
}