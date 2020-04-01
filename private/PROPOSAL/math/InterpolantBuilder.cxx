
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/math/Interpolant.h"
#include <memory>
using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Defaults for InterpolantBuilder
// ------------------------------------------------------------------------- //

const int InterpolantBuilder::default_max = 1.0;
const double InterpolantBuilder::default_xmin = 1.0;
const double InterpolantBuilder::default_xmax = 1.0;

const int InterpolantBuilder::default_romberg = 1.0;
const bool InterpolantBuilder::default_rational = false;
const bool InterpolantBuilder::default_relative = false;
const bool InterpolantBuilder::default_isLog = false;

const int InterpolantBuilder::default_rombergY = 1.0;
const bool InterpolantBuilder::default_rationalY = false;
const bool InterpolantBuilder::default_relativeY = false;
const bool InterpolantBuilder::default_logSubst = false;

const Interpolant2DBuilder::Function2D Interpolant2DBuilder::default_function2d
    = NULL;

// ------------------------------------------------------------------------- //
// Interpolant1D Builder
// ------------------------------------------------------------------------- //

Interpolant1DBuilder::Definition::Definition(Function1D func, int max, double xmin,
    double xmax, int romberg, bool rational, bool relative, bool isLog,
    int rombergY, bool rationalY, bool relativeY, bool logSubst)
    : function1d(func)
    , max(max)
    , xmin(xmin)
    , xmax(xmax)
    , romberg(romberg)
    , rational(rational)
    , relative(relative)
    , isLog(isLog)
    , rombergY(rombergY)
    , rationalY(rationalY)
    , relativeY(relativeY)
    , logSubst(logSubst)
{
}

Interpolant1DBuilder::Interpolant1DBuilder()
    : builder_def(new Definition())
{
}

Interpolant1DBuilder::Interpolant1DBuilder(const Definition& builder_def)
    : builder_def(new Interpolant1DBuilder::Definition(builder_def))
{
}

Interpolant1DBuilder::~Interpolant1DBuilder() = default;

Interpolant1DBuilder& Interpolant1DBuilder::SetFunction1D(Function1D val)
{
    builder_def->function1d = val;
    return *this;
}

Interpolant1DBuilder& Interpolant1DBuilder::SetMax(const int val)
{
    builder_def->max = val;
    return *this;
}
Interpolant1DBuilder& Interpolant1DBuilder::SetXMin(const double val)
{
    builder_def->xmin = val;
    return *this;
}
Interpolant1DBuilder& Interpolant1DBuilder::SetXMax(const double val)
{
    builder_def->xmax = val;
    return *this;
}

Interpolant1DBuilder& Interpolant1DBuilder::SetRomberg(const int val)
{
    builder_def->romberg = val;
    return *this;
}
Interpolant1DBuilder& Interpolant1DBuilder::SetRational(const bool val)
{
    builder_def->rational = val;
    return *this;
}
Interpolant1DBuilder& Interpolant1DBuilder::SetRelative(const bool val)
{
    builder_def->relative = val;
    return *this;
}
Interpolant1DBuilder& Interpolant1DBuilder::SetIsLog(const bool val)
{
    builder_def->isLog = val;
    return *this;
}

Interpolant1DBuilder& Interpolant1DBuilder::SetRombergY(const int val)
{
    builder_def->rombergY = val;
    return *this;
}
Interpolant1DBuilder& Interpolant1DBuilder::SetRationalY(const bool val)
{
    builder_def->rationalY = val;
    return *this;
}
Interpolant1DBuilder& Interpolant1DBuilder::SetRelativeY(const bool val)
{
    builder_def->relativeY = val;
    return *this;
}
Interpolant1DBuilder& Interpolant1DBuilder::SetLogSubst(const bool val)
{
    builder_def->logSubst = val;
    return *this;
}

std::unique_ptr<Interpolant> Interpolant1DBuilder::build()
{
    std::unique_ptr<Interpolant> interpolant(new Interpolant(builder_def->max,
        builder_def->xmin, builder_def->xmax, builder_def->function1d,
        builder_def->romberg, builder_def->rational, builder_def->relative,
        builder_def->isLog, builder_def->rombergY, builder_def->rationalY,
        builder_def->relativeY, builder_def->logSubst));
    return interpolant;
}

// ------------------------------------------------------------------------- //
// Interpolant2D Builder
// ------------------------------------------------------------------------- //

Interpolant2DBuilder::Interpolant2DBuilder()
    : InterpolantBuilder()
    , function2d(default_function2d)
    , max1(default_max)
    , x1min(default_xmin)
    , x1max(default_xmax)
    , max2(default_max)
    , x2min(default_xmin)
    , x2max(default_xmax)
    , romberg1(default_romberg)
    , rational1(default_rational)
    , relative1(default_relative)
    , isLog1(default_isLog)
    , romberg2(default_romberg)
    , rational2(default_rational)
    , relative2(default_relative)
    , isLog2(default_isLog)
    , rombergY(default_rombergY)
    , rationalY(default_rationalY)
    , relativeY(default_relativeY)
    , logSubst(default_logSubst)
{
}

Interpolant2DBuilder::Interpolant2DBuilder(const Interpolant2DBuilder& builder)
    : function2d(builder.function2d)
    , max1(builder.max1)
    , x1min(builder.x1min)
    , x1max(builder.x1max)
    , max2(builder.max2)
    , x2min(builder.x2min)
    , x2max(builder.x2max)
    , romberg1(builder.romberg1)
    , rational1(builder.rational1)
    , relative1(builder.relative1)
    , isLog1(builder.isLog1)
    , romberg2(builder.romberg2)
    , rational2(builder.rational2)
    , relative2(builder.relative2)
    , isLog2(builder.isLog2)
    , rombergY(builder.rombergY)
    , rationalY(builder.rationalY)
    , relativeY(builder.relativeY)
    , logSubst(builder.logSubst)
{
}

Interpolant2DBuilder_array_as::Interpolant2DBuilder_array_as()
    : InterpolantBuilder()
    , romberg1(default_romberg)
    , rational1(default_rational)
    , relative1(default_relative)
    , romberg2(default_romberg)
    , rational2(default_rational)
    , relative2(default_relative)
{
}

Interpolant2DBuilder_array_as::Interpolant2DBuilder_array_as(
    const Interpolant2DBuilder_array_as& builder)
    : x1(builder.x1)
    , x2(builder.x2)
    , y(builder.y)
    , romberg1(builder.romberg1)
    , rational1(builder.rational1)
    , relative1(builder.relative1)
    , romberg2(builder.romberg2)
    , rational2(builder.rational2)
    , relative2(builder.relative2)
{
}

std::unique_ptr<Interpolant> Interpolant2DBuilder::build()
{
    return std::unique_ptr<Interpolant>(
        new Interpolant(max1, x1min, x1max, max2, x2min, x2max, function2d,
            romberg1, rational1, relative1, isLog1, romberg2, rational2,
            relative2, isLog2, rombergY, rationalY, relativeY, logSubst));
}

std::unique_ptr<Interpolant> Interpolant2DBuilder_array_as::build()
{
    return std::unique_ptr<Interpolant>(new Interpolant(x1, x2, y, romberg1,
        rational1, relative1, romberg2, rational2, relative2));
}
