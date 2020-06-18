#pragma once

#include "PROPOSAL/crossection/CrossSectionDNDX/CrossSectionDNDX.h"
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

    double Calculate(
        double energy, double v, v_trafo_t trafo = nullptr) override
    {
        auto integral_lim = GetIntegrationLimits(energy);
        if (trafo)
            v = trafo(get<MIN>(integral_lim), get<MAX>(integral_lim), v);
        return dndx_integral(integral, energy, get<MIN>(integral_lim), v, 0);
    }

    double GetUpperLimit(
        double energy, double rate, v_trafo_t trafo = nullptr) override
    {
        auto integral_lim = GetIntegrationLimits(energy);
        dndx_integral(integral, energy, get<MIN>(integral_lim),
            get<MAX>(integral_lim), -rate);
        auto v = integral.GetUpperLimit();
        if (trafo)
            v = trafo(get<MIN>(integral_lim), get<MAX>(integral_lim), v);
        return v;
    }
};
} // namespace PROPOSAL
