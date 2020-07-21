#pragma once

#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXIntegral.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include <type_traits>

using std::unique_ptr;

namespace PROPOSAL {

double transform_relativ_loss(double v_cut, double v_max, double v);
double retransform_relativ_loss(double v_cut, double v_max, double v);

class CrossSectionDNDXInterpolant : public CrossSectionDNDXIntegral {
    unique_ptr<Interpolant> dndx;

    unique_ptr<Interpolant> build_dndx(
        Interpolant2DBuilder::Definition, const InterpolationDef&);

public:
    template <typename Param, typename Particle, typename Target>
    CrossSectionDNDXInterpolant(const Param& _param, const Particle& _p_def,
        const Target& _target, shared_ptr<const EnergyCutSettings> _cut,
        const InterpolationDef& _interpol_def)
        : CrossSectionDNDXIntegral(_param, _p_def, _target, _cut)
        , dndx(
              build_dndx(build_dndx_interpol_def(_param, _p_def, _interpol_def),
                  _interpol_def))
    {
    }

    double Calculate(double) final;
    double Calculate(double, double) final;
    double GetUpperLimit(double, double) final;
};

Interpolant2DBuilder::Definition build_dndx_interpol_def(
    const crosssection::Parametrization&, const ParticleDef&, const InterpolationDef&);
} // namespace PROPOSAL
