
#include <boost/bind.hpp>

#include "PROPOSAL/crossection/BremsInterpolant.h"
#include "PROPOSAL/crossection/BremsIntegral.h"
#include "PROPOSAL/crossection/parametrization/Parametrization.h"

#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/methods.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Output.h"

using namespace PROPOSAL;

BremsInterpolant::BremsInterpolant(Parametrization& param): CrossSectionInterpolant(param)
{
    Parametrization::Definition param_def = parametrization_.GetDefinition();
    Interpolant1DBuilder builder1d;
    Helper::InterpolantBuilderContainer builder_container;

    // Needed for CalculatedEdx integration
    BremsIntegral brems(param);

    builder1d.SetMax(NUM1)
        .SetXMin(param.GetParticleDef().low)
        .SetXMax(BIGENERGY)
        .SetRomberg(param_def.order_of_interpolation)
        .SetRational(true)
        .SetRelative(false)
        .SetIsLog(true)
        .SetRombergY(param_def.order_of_interpolation)
        .SetRationalY(false)
        .SetRelativeY(false)
        .SetLogSubst(true)
        .SetFunction1D(boost::bind(&CrossSection::CalculatedEdx, brems, _1));

    builder_container.push_back(std::make_pair(&builder1d, &dedx_interpolant_));

    log_info("Initialize dEdx for %s", typeid(parametrization_).name());
    Helper::InitializeInterpolation("dEdx",
                                    builder_container,
                                    std::vector<Parametrization*>(1, &parametrization_));
    log_info("Initialization dEdx for %s done!", typeid(parametrization_).name());
}

BremsInterpolant::BremsInterpolant(const BremsInterpolant& brems): CrossSectionInterpolant(brems)
{
}

BremsInterpolant::~BremsInterpolant()
{
}

// ----------------------------------------------------------------- //
// Public methods
// ----------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double BremsInterpolant::CalculatedEdx(double energy)
{
    if (parametrization_.GetMultiplier() <= 0)
    {
        return 0;
    }

    return std::max(dedx_interpolant_->Interpolate(energy), 0.0);
}
