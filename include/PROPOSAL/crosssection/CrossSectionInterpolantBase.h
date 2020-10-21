#pragma once

#include "PROPOSAL/math/InterpolantBuilder.h"

#include <memory>

namespace PROPOSAL {

struct CrossSectionInterpolantBase {
    static std::unique_ptr<Interpolant2DBuilder::Definition> dNdx_def;
    static std::unique_ptr<Interpolant1DBuilder::Definition> dEdx_def;
    static std::unique_ptr<Interpolant1DBuilder::Definition> dE2dx_def;

    Interpolant2DBuilder::Definition const& GetDNdxDef();
    Interpolant1DBuilder::Definition const& GetDEdxDef();
    Interpolant1DBuilder::Definition const& GetDE2dxDef();
    void SetDNdxDef(Interpolant2DBuilder::Definition const& def);
    void SetDEdxDef(Interpolant1DBuilder::Definition const& def);
    void SetDE2dxDef(Interpolant1DBuilder::Definition const& def);
};


} // namespace PROPOSAL
