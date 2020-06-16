#pragma once

#include "PROPOSAL/crossection/CrossSectionDNDX/CrossSectionDNDXIntegral.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/crossection/CrossSectionInterpolant.h"

#include <type_traits>

using std::unique_ptr;

namespace PROPOSAL {

class CrossSectionDNDXInterpolant : public CrossSectionDNDXIntegral {
    unique_ptr<Interpolant> dndx;

    unique_ptr<Interpolant> build_dndx(
        Interpolant2DBuilder::Definition, const InterpolationDef&);

public:
    template <typename Param, typename Particle, typename Target>
    CrossSectionDNDXInterpolant(const Param& _param, const Particle& _p_def,
        const Target& _target, shared_ptr<EnergyCutSettings> _cut,
        const InterpolationDef& _interpol_def)
        : CrossSectionDNDXIntegral(_param, _p_def, _target, _cut)
        , dndx(
              build_dndx(build_dndx_interpol_def(_param, _p_def, _interpol_def),
                  _interpol_def))
    {
    }

    double Calculate(double, double, v_trafo_t = nullptr) final;
    double GetUpperLim(double, double, v_trafo_t = nullptr) final;
};

Interpolant2DBuilder::Definition build_dndx_interpol_def(
    const Parametrization&, const ParticleDef&, const InterpolationDef&);
} // namespace PROPOSAL
