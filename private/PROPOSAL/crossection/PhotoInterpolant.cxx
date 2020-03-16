
#include <functional>

#include "PROPOSAL/crossection/PhotoIntegral.h"
#include "PROPOSAL/crossection/PhotoInterpolant.h"
#include "PROPOSAL/crossection/parametrization/Photonuclear.h"

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

PhotoInterpolant::PhotoInterpolant(const Photonuclear& param, std::shared_ptr<const EnergyCutSettings> cuts, const InterpolationDef& def)
    : CrossSectionInterpolant(param, cuts)
{
    // Use parent CrossSection dNdx interpolation
    InitdNdxInterpolation(def);

    // --------------------------------------------------------------------- //
    // Builder for DEdx
    // --------------------------------------------------------------------- //

    Interpolant1DBuilder builder1d;
    Helper::InterpolantBuilderContainer builder_container;

    // Needed for CalculatedEdx integration
    PhotoIntegral photo(param, cuts);

    builder1d.SetMax(def.nodes_cross_section)
        .SetXMin(param.GetParticleMass())
        .SetXMax(def.max_node_energy)
        .SetRomberg(def.order_of_interpolation)
        .SetRational(true)
        .SetRelative(false)
        .SetIsLog(true)
        .SetRombergY(def.order_of_interpolation)
        .SetRationalY(false)
        .SetRelativeY(false)
        .SetLogSubst(false)
        .SetFunction1D(std::bind(&CrossSectionIntegral::CalculatedEdxWithoutMultiplier, &photo, std::placeholders::_1));

    builder_container.push_back(&builder1d);

    // --------------------------------------------------------------------- //
    // Builder for DE2dx
    // --------------------------------------------------------------------- //

    Interpolant1DBuilder builder_de2dx;
    Helper::InterpolantBuilderContainer builder_container_de2dx;

    builder_de2dx.SetMax(def.nodes_continous_randomization)
        .SetXMin(param.GetParticleMass())
        .SetXMax(def.max_node_energy)
        .SetRomberg(def.order_of_interpolation)
        .SetRational(false)
        .SetRelative(false)
        .SetIsLog(true)
        .SetRombergY(def.order_of_interpolation)
        .SetRationalY(false)
        .SetRelativeY(false)
        .SetLogSubst(false)
        .SetFunction1D(std::bind(&PhotoIntegral::CalculatedE2dxWithoutMultiplier, &photo, std::placeholders::_1));

    builder_container_de2dx.push_back(&builder_de2dx);

    auto dedx_interpolant_vec = Helper::InitializeInterpolation("dEdx", builder_container, std::vector<Parametrization*>(1, parametrization_), def);
    dedx_interpolant_ = std::move(dedx_interpolant_vec.at(0));

    auto de2dx_interpolant_vec = Helper::InitializeInterpolation(
        "dE2dx", builder_container_de2dx, std::vector<Parametrization*>(1, parametrization_), def);
    de2dx_interpolant_ = std::move(de2dx_interpolant_vec.at(0));
}

/*PhotoInterpolant::PhotoInterpolant(const PhotoInterpolant& photo)
    : CrossSectionInterpolant(photo)
{
}*/

PhotoInterpolant::~PhotoInterpolant() {}

// ----------------------------------------------------------------- //
// Public methods
// ----------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double PhotoInterpolant::CalculatedEdx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    return parametrization_->GetMultiplier() * std::max(dedx_interpolant_->Interpolate(energy), 0.0);
}
