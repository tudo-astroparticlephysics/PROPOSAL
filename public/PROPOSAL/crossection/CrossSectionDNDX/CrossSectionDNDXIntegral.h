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
    CrossSectionDNDXIntegral(Param _param, Particle _p_def, Target _target,
        shared_ptr<EnergyCutSettings> cut)
        : CrossSectionDNDX(_param, _p_def, _target, cut)
    {
    }

    double Calculate(double, double, v_trafo_t = nullptr) override;
    double GetUpperLim(double, double, v_trafo_t = nullptr) override;
};
} // namespace PROPOSAL
