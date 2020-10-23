#pragma once

#include "PROPOSAL/crosssection/CrossSectionInterpolantBase.h"
#include "PROPOSAL/crosssection/CrossSectionDNDX/CrossSectionDNDXIntegral.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include <type_traits>

using std::unique_ptr;

namespace PROPOSAL {
namespace crosssection{
    struct Parametrization;
}
struct ParticleDef;

double transform_relativ_loss(double v_cut, double v_max, double v);
double retransform_relativ_loss(double v_cut, double v_max, double v);

class CrossSectionDNDXInterpolant : public CrossSectionDNDXIntegral , public CrossSectionInterpolantBase {
    unique_ptr<Interpolant> dndx;
    unique_ptr<Interpolant> dndx1d;

    unique_ptr<Interpolant> build_dndx(crosssection::Parametrization const&, ParticleDef const&);
    unique_ptr<Interpolant> build_dndx1d(crosssection::Parametrization const&, ParticleDef const&);
public:
    template <typename Param, typename Particle, typename Target>
    CrossSectionDNDXInterpolant(const Param& _param, const Particle& _p_def, const Target& _target, shared_ptr<const EnergyCutSettings> _cut)
        : CrossSectionDNDXIntegral(_param, _p_def, _target, _cut)
        , dndx(build_dndx(_param, _p_def))
        , dndx1d(build_dndx1d(_param, _p_def))
    {
    }


    double Calculate(double) final;
    double Calculate(double, double) final;
    double GetUpperLimit(double, double) final;
};
} // namespace PROPOSAL
