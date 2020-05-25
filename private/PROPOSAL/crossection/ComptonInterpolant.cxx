
#include <cmath>
#include <functional>

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

/* ComptonInterpolant::ComptonInterpolant(unique_ptr<Compton>&& param, */
/*     shared_ptr<const EnergyCutSettings> cut, const InterpolationDef& def) */
/*     : CrossSectionInterpolant(forward<unique_ptr<Compton>>(param), cut, def) */
/* { */
/* } */

/* unordered_map<size_t, unique_ptr<Interpolant>> ComptonInterpolant::init_dndx_interpolation( */
/*     const InterpolationDef& def) */
/* { */
/*     Interpolant2DBuilder::Definition interpolant_def; */
/*     interpolant_def.max1 = def.nodes_cross_section; */
/*     interpolant_def.x1min = ME; */
/*     interpolant_def.x1max = def.max_node_energy; */
/*     interpolant_def.max2 = def.nodes_cross_section; */
/*     interpolant_def.x2min = 1. / (2. * (1. - def.nodes_cross_section)); */
/*     interpolant_def.x2max = (1. - 2. * def.nodes_cross_section) */
/*         / (2. * (1 - def.nodes_cross_section)); */
/*     interpolant_def.romberg1 = def.order_of_interpolation; */
/*     interpolant_def.isLog1 = true; */
/*     interpolant_def.romberg2 = def.order_of_interpolation; */
/*     interpolant_def.rational2 = true; */
/*     interpolant_def.rombergY = def.order_of_interpolation; */
/*     interpolant_def.rationalY = true; */

/*     unordered_map<size_t, unique_ptr<Interpolant>> dndx_; */
/*     for (auto dndx : dndx_integral_) { */
/*         interpolant_def.function2d = dndx.second; */

/*         auto builder = make_unique<Interpolant2DBuilder>(interpolant_def); */
/*         auto interpolant = Helper::InitializeInterpolation( */
/*             "dNdx", std::move(builder), GetHash(), def); */

/*         dndx_[dndx.first] = std::move(interpolant); */
/*     } */

/*     return dndx_; */
/* } */
