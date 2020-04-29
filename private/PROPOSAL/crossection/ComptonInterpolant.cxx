
#include <cmath>
#include <functional>

#include "PROPOSAL/crossection/ComptonIntegral.h"
#include "PROPOSAL/crossection/ComptonInterpolant.h"
#include "PROPOSAL/crossection/parametrization/Compton.h"

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;
using std::placeholders::_1;

vector<unique_ptr<Interpolant>> ComptonInterpolant::init_dndx_interpolation(
    const InterpolationDef& def)
{
    Interpolant2DBuilder::Definition interpolant_def;
    interpolant_def.max1 = def.nodes_cross_section;
    interpolant_def.x1min = ME;
    interpolant_def.x1max = def.max_node_energy;
    interpolant_def.max2 = def.nodes_cross_section;
    interpolant_def.x2min = 1. / (2. * (1. - def.nodes_cross_section));
    interpolant_def.x2max = (1. - 2. * def.nodes_cross_section)
        / (2. * (1 - def.nodes_cross_section));
    interpolant_def.romberg1 = def.order_of_interpolation;
    interpolant_def.isLog1 = true;
    interpolant_def.romberg2 = def.order_of_interpolation;
    interpolant_def.rational2 = true;
    interpolant_def.rombergY = def.order_of_interpolation;
    interpolant_def.rationalY = true;

    vector<unique_ptr<InterpolantBuilder>> builders;
    for (auto dndx : dndx_integral_) {
        interpolant_def.function2d = dndx;
        builders.emplace_back(new Interpolant2DBuilder(interpolant_def));
    }

    const auto& hash = parametrization_->GetHash();
    auto aux = Helper::InitializeInterpolation("dNdx_diff", builders, hash, def);
    return aux;
}
