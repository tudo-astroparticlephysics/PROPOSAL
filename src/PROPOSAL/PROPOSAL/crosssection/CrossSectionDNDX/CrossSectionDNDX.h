#pragma once

#include "PROPOSAL/crosssection/parametrization/Parametrization.h"

#include <functional>
#include <spdlog/fwd.h>

namespace PROPOSAL {
class EnergyCutSettings;
}

namespace PROPOSAL {
class CrossSectionDNDX {
    using param_medium_t = crosssection::Parametrization<Medium>;
    using param_comp_t = crosssection::Parametrization<Component>;
    using lim_func_t = std::function<crosssection::KinematicLimits(double)>;

protected:
    size_t hash;
    std::shared_ptr<spdlog::logger> logger;
    double lower_energy_lim;

private:
    lim_func_t kinematic_limits;
    std::shared_ptr<const EnergyCutSettings> cut;

    CrossSectionDNDX(lim_func_t kin_lim, double lower_energy_lim,
        std::shared_ptr<const EnergyCutSettings> cut, size_t hash);

public:
    CrossSectionDNDX(param_medium_t const&, ParticleDef, Medium,
        std::shared_ptr<const EnergyCutSettings>, size_t hash);

    CrossSectionDNDX(param_comp_t const&, ParticleDef, Component,
        std::shared_ptr<const EnergyCutSettings>, size_t hash);

    virtual ~CrossSectionDNDX() = default;

    virtual double Calculate(double energy) = 0;
    virtual double Calculate(double energy, double v) = 0;
    virtual double GetUpperLimit(double energy, double rate) = 0;

    struct IntegrationLimit {
        double min, max;
    };
    IntegrationLimit GetIntegrationLimits(double energy) const;
    double GetLowerEnergyLim() const;
    size_t GetHash() const noexcept { return hash; }
};
} // namespace PROPOSAL
