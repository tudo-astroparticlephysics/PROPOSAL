
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/math/Interpolant.h"

#include <memory>

using std::unique_ptr;
using namespace PROPOSAL;

unique_ptr<Interpolant> Interpolant1DBuilder::build()
{
    return unique_ptr<Interpolant>(new Interpolant(builder_def->max,
        builder_def->xmin, builder_def->xmax, builder_def->function1d,
        builder_def->romberg, builder_def->rational, builder_def->relative,
        builder_def->isLog, builder_def->rombergY, builder_def->rationalY,
        builder_def->relativeY, builder_def->logSubst));
}

unique_ptr<Interpolant> Interpolant2DBuilder::build()
{
    return unique_ptr<Interpolant>(new Interpolant(builder_def->max1,
        builder_def->x1min, builder_def->x1max, builder_def->max2,
        builder_def->x2min, builder_def->x2max, builder_def->function2d,
        builder_def->romberg1, builder_def->rational1, builder_def->relative1,
        builder_def->isLog1, builder_def->romberg2, builder_def->rational2,
        builder_def->relative2, builder_def->isLog2, builder_def->rombergY,
        builder_def->rationalY, builder_def->relativeY, builder_def->logSubst));
}
