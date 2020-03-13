
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/math/Interpolant.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Defaults for InterpolantBuilder
// ------------------------------------------------------------------------- //

const int InterpolantBuilder::default_max     = 1.0;
const double InterpolantBuilder::default_xmin = 1.0;
const double InterpolantBuilder::default_xmax = 1.0;

const int InterpolantBuilder::default_romberg   = 1.0;
const bool InterpolantBuilder::default_rational = false;
const bool InterpolantBuilder::default_relative = false;
const bool InterpolantBuilder::default_isLog    = false;

const int InterpolantBuilder::default_rombergY   = 1.0;
const bool InterpolantBuilder::default_rationalY = false;
const bool InterpolantBuilder::default_relativeY = false;
const bool InterpolantBuilder::default_logSubst  = false;

const Interpolant1DBuilder::Function1D Interpolant1DBuilder::default_function1d = NULL;
const Interpolant2DBuilder::Function2D Interpolant2DBuilder::default_function2d = NULL;

// ------------------------------------------------------------------------- //
// Interpolant1D Builder
// ------------------------------------------------------------------------- //

Interpolant1DBuilder::Interpolant1DBuilder()
    : InterpolantBuilder()
    , function1d(default_function1d)
    , max(default_max)
    , xmin(default_xmin)
    , xmax(default_xmax)
    , romberg(default_romberg)
    , rational(default_rational)
    , relative(default_relative)
    , isLog(default_isLog)
    , rombergY(default_rombergY)
    , rationalY(default_rationalY)
    , relativeY(default_relativeY)
    , logSubst(default_logSubst)
{
}

Interpolant1DBuilder::Interpolant1DBuilder(const Interpolant1DBuilder& builder)
    : function1d(builder.function1d)
    , max(builder.max)
    , xmin(builder.xmin)
    , xmax(builder.xmax)
    , romberg(builder.romberg)
    , rational(builder.rational)
    , relative(builder.relative)
    , isLog(builder.isLog)
    , rombergY(builder.rombergY)
    , rationalY(builder.rationalY)
    , relativeY(builder.relativeY)
    , logSubst(builder.logSubst)
{
}

Interpolant* Interpolant1DBuilder::build()
{
    return new Interpolant(
        max, xmin, xmax, function1d, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, logSubst);
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

Interpolant2DBuilder_array_as::Interpolant2DBuilder_array_as(const Interpolant2DBuilder_array_as& builder)
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

Interpolant* Interpolant2DBuilder::build()
{
    return new Interpolant(max1,
                           x1min,
                           x1max,
                           max2,
                           x2min,
                           x2max,
                           function2d,
                           romberg1,
                           rational1,
                           relative1,
                           isLog1,
                           romberg2,
                           rational2,
                           relative2,
                           isLog2,
                           rombergY,
                           rationalY,
                           relativeY,
                           logSubst);
}

Interpolant* Interpolant2DBuilder_array_as::build()
{
    return new Interpolant(x1,
            x2,
            y,
            romberg1,
            rational1,
            relative1,
            romberg2,
            rational2,
            relative2);
}

