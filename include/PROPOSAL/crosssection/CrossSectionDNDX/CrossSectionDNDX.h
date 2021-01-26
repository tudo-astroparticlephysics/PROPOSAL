#pragma once

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/crosssection/parametrization/Parametrization.h"

#include <memory>

namespace PROPOSAL {
class CrossSectionDNDX {
    std::function<std::tuple<double, double>(double)> kinematic_limits;
    std::shared_ptr<const EnergyCutSettings> cut;
    size_t hash;

    template <typename T1, typename T2, typename T3>
    inline auto define_kinematic_limits(T1 param, T2 particle, T3 target)
    {
        return [param, particle, target](double E) {
            return param.GetKinematicLimits(particle, target, E);
        };
    }

public:
    template <typename T1, typename T2, typename T3>
    CrossSectionDNDX(T1 _param, T2 _particle, T3 _target)
        : kinematic_limits(define_kinematic_limits(_param, _particle, _target))
        , hash(0)
    {
        auto logger = Logging::Get("PROPOSAL.CrossSection.dNdx");
        logger->info("Building {} dNdx.", _param.name);
        hash_combine(
            hash, _param.GetHash(), _particle.GetHash(), _target.GetHash());
    }

    template <typename T1, typename T2, typename T3>
    CrossSectionDNDX(T1&& _param, T2&& _particle, T3&& _target,
        EnergyCutSettings const& _cut)
        : CrossSectionDNDX(std::forward<T1>(_param),
            std::forward<T2>(_particle), std::forward<T3>(_target))
    {
        cut = std::make_shared<const EnergyCutSettings>(_cut);
        hash_combine(hash, _cut.GetHash());
    }

    virtual ~CrossSectionDNDX() = default;

    virtual double Calculate(double energy) = 0;
    virtual double Calculate(double energy, double v) = 0;
    virtual double GetUpperLimit(double energy, double rate) = 0;

    enum { MIN, MAX };
    std::array<double, 2> GetIntegrationLimits(double energy) const;

    size_t GetHash() const noexcept { return hash; }
};
} // namespace PROPOSAL
