
#pragma once

#include "PROPOSAL/crossection/CrossSectionIntegral.h"

namespace PROPOSAL {

class Bremsstrahlung;

class BremsIntegral : public CrossSectionIntegral
{
public:
    BremsIntegral(const Bremsstrahlung&);
    BremsIntegral(const BremsIntegral&);
    virtual ~BremsIntegral();

    CrossSection* clone() const { return new BremsIntegral(*this); }

    // ----------------------------------------------------------------- //
    // Public methods
    // ----------------------------------------------------------------- //

    double CalculatedEdx(double energy);
};

} // namespace PROPOSAL
