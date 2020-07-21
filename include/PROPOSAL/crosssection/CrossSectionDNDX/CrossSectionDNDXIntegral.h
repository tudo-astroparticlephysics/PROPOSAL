#pragma once

#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDX.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"

#include <memory>
#include <unordered_map>

using std::shared_ptr;
using std::unordered_map;

using PROPOSAL::Components::Component;

namespace PROPOSAL {
class CrossSectionDNDXIntegral : public CrossSectionDNDX {
    Integral integral;

public:
    template <typename Param, typename Particle, typename Target>
    CrossSectionDNDXIntegral(Param _param, Particle _particle, Target _target,
        shared_ptr<const EnergyCutSettings> _cut)
        : CrossSectionDNDX(_param, _particle, _target, _cut)
    {
    }

    double Calculate(double energy) override
    {
        return Calculate(energy, get<MAX>(GetIntegrationLimits(energy)));
    }

    double Calculate(double energy, double v) override
    {
        auto integral_lim = GetIntegrationLimits(energy);
        return dndx_integral(integral, energy, get<MIN>(integral_lim), v, 0);
    }

    double GetUpperLimit(double energy, double rate) override
    {
        auto integral_lim = GetIntegrationLimits(energy);
        dndx_integral(integral, energy, get<MIN>(integral_lim),
            get<MAX>(integral_lim), -rate);
        auto v = integral.GetUpperLimit();
        return v;
    }
};
} // namespace PROPOSAL
