#include "PROPOSAL/crosssection/CrossSectionInterpolantBase.h"

namespace PROPOSAL {
    std::unique_ptr<Interpolant2DBuilder::Definition> CrossSectionInterpolantBase::dNdx_def = nullptr;
    std::unique_ptr<Interpolant1DBuilder::Definition> CrossSectionInterpolantBase::dEdx_def = nullptr;
    std::unique_ptr<Interpolant1DBuilder::Definition> CrossSectionInterpolantBase::dE2dx_def = nullptr;


    Interpolant2DBuilder::Definition const& CrossSectionInterpolantBase::GetDNdxDef() {
        if(not dNdx_def)
            dNdx_def = std::make_unique<Interpolant2DBuilder::Definition>();
        return *dNdx_def;
    }

    Interpolant1DBuilder::Definition const& CrossSectionInterpolantBase::GetDEdxDef() {
        if(not dEdx_def)
            dEdx_def = std::make_unique<Interpolant1DBuilder::Definition>();
        return *dEdx_def;
    }

    Interpolant1DBuilder::Definition const& CrossSectionInterpolantBase::GetDE2dxDef() {
        if(not dE2dx_def)
            dE2dx_def = std::make_unique<Interpolant1DBuilder::Definition>();
        return *dE2dx_def;
    }

    void CrossSectionInterpolantBase::SetDNdxDef(Interpolant2DBuilder::Definition const& def) {
        dNdx_def = std::make_unique<Interpolant2DBuilder::Definition>(def);
    }

    void CrossSectionInterpolantBase::SetDEdxDef(Interpolant1DBuilder::Definition const& def) {
        dEdx_def = std::make_unique<Interpolant1DBuilder::Definition>(def);
    }

    void CrossSectionInterpolantBase::SetDE2dxDef(Interpolant1DBuilder::Definition const& def) {
        dE2dx_def = std::make_unique<Interpolant1DBuilder::Definition>(def);
    }
} // namespace PROPOSAL
