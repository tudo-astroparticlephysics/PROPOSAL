#pragma once

#include <utility>
#include <nlohmann/json.hpp>
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/scattering/stochastic_deflection/Parametrization.h"
#include "PROPOSAL/scattering/multiple_scattering/Parametrization.h"
#include "PROPOSAL/scattering/stochastic_deflection/ScatteringFactory.h"
#include "PROPOSAL/scattering/multiple_scattering/ScatteringFactory.h"
#include "PROPOSAL/scattering/ScatteringMultiplier.h"

namespace PROPOSAL {

template <typename Cross>
std::unique_ptr<Scattering> make_scattering(const nlohmann::json& config,
    ParticleDef const& p_def,Medium const& medium, Cross&& cross, bool interpol) {
    std::unique_ptr<multiple_scattering::Parametrization> ms_pointer = nullptr;
    if (config.contains("multiple_scattering"))
        ms_pointer = make_multiple_scattering(config["multiple_scattering"].get<std::string>(),
                                              p_def, medium, cross, interpol);
    double ms_multiplier = config.value("multiple_scattering_multiplier", 1.0);

    std::vector<std::unique_ptr<stochastic_deflection::Parametrization>> sd = {};
    if (config.contains("stochastic_deflection")) {
        assert(config["stochastic_deflection"].is_array());
        for (const auto& param_raw : config.at("stochastic_deflection")) {
            std::string param = param_raw;
            std::transform(param.begin(), param.end(), param.begin(), ::tolower);
            try {
                // if param matches the name of a parametrization, use it
                sd.emplace_back(make_stochastic_deflection(param, p_def, medium));
            } catch (const std::out_of_range&) {
                if (param.rfind("default", 0) == 0) {
                    // if param starts with "default", find default stochastic
                    // deflection corresponding to interaction type
                    param.erase(0,7);
                    auto it = std::find_if(std::begin(Type_Interaction_Name_Map),
                                           std::end(Type_Interaction_Name_Map),
                                           [& param](auto&& p) { std::string type_lower = p.second;
                                           std::transform(type_lower.begin(), type_lower.end(), type_lower.begin(), ::tolower);
                                           return type_lower == param;}
                                           );
                    if (it == std::end(Type_Interaction_Name_Map)) {
                        std::ostringstream s;
                        s << "Interaction type " << param << " not found.";
                        throw std::out_of_range(s.str());
                    } else {
                        sd.emplace_back(make_default_stochastic_deflection(it->first, p_def, medium));
                    }
                } else {
                    throw;
                }
            }
        }
    }

    std::vector<std::pair<InteractionType, double>> sd_multiplier = {};
    if (config.contains("stochastic_deflection_multiplier")) {
        assert(config["stochastic_deflection_multiplier"].is_object());
        for (auto& item : config["stochastic_deflection_multiplier"].items()) {
            std::string type = item.key();
            std::transform(type.begin(), type.end(), type.begin(), ::tolower);
            auto it = std::find_if(std::begin(Type_Interaction_Name_Map),
                                   std::end(Type_Interaction_Name_Map),
                                   [& type](auto&& p) { std::string type_lower = p.second;
                                       std::transform(type_lower.begin(), type_lower.end(), type_lower.begin(), ::tolower);
                                       return type_lower == type;});
            if (it == std::end(Type_Interaction_Name_Map)) {
                std::ostringstream s;
                s << "Interaction type " << type << " not found.";
                throw std::out_of_range(s.str());
            }
            sd_multiplier.push_back(std::make_pair(it->first, item.value()));
        }
    }

    if (ms_multiplier == 1.0 && sd_multiplier.empty()) {
        return std::make_unique<Scattering>(ms_pointer, sd);
    } else {
        return std::unique_ptr<Scattering>(new ScatteringMultiplier(ms_pointer, sd, ms_multiplier, sd_multiplier));
    }

}
}
