
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/methods.h"

#include <memory>
#include <iostream>

using namespace PROPOSAL;

size_t Interpolant1DBuilder::Definition::GetHash() {
    auto hash = (size_t)0;
    hash_combine(hash, nodes, xmin, xmax, romberg,
            rational, relative, isLog, rombergY, rationalY,
            relativeY, logSubst);
    return hash;
}

Interpolant1DBuilder::Interpolant1DBuilder(Interpolant1DBuilder::Definition _def)
    : def(std::make_unique<Interpolant1DBuilder::Definition>(_def))
{
}

std::unique_ptr<Interpolant> Interpolant1DBuilder::build() const
{
    return std::make_unique<Interpolant>(def->nodes, def->xmin, def->xmax,
            def->function1d, def->romberg, def->rational, def->relative, def->isLog,
            def->rombergY, def->rationalY, def->relativeY, def->logSubst);
}

size_t Interpolant2DBuilder::Definition::GetHash() {
    auto hash = Interpolant1DBuilder::Definition::GetHash();
    hash_combine(hash, nodes2, x2min, x2max, romberg2,
        rational2, relative2, isLog2);
    return hash;
}

Interpolant2DBuilder::Interpolant2DBuilder(Interpolant2DBuilder::Definition _def)
    : def(std::make_unique<Interpolant2DBuilder::Definition>(_def))
{
}

std::unique_ptr<Interpolant> Interpolant2DBuilder::build() const
{
    return std::make_unique<Interpolant>(def->nodes,
            def->xmin, def->xmax, def->nodes2,
            def->x2min, def->x2max, def->function2d,
            def->romberg, def->rational, def->relative,
            def->isLog, def->romberg2, def->rational2,
            def->relative2, def->isLog2, def->rombergY,
            def->rationalY, def->relativeY, def->logSubst);
}
