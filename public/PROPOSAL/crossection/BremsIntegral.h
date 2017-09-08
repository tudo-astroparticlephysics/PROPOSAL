
#pragma once


#include "PROPOSAL/crossection/CrossSectionIntegral.h"

namespace PROPOSAL
{

class BremsIntegral: public CrossSectionIntegral
{
    public:
        BremsIntegral(Parametrization&);
        BremsIntegral(const BremsIntegral&);
        virtual ~BremsIntegral();

        CrossSection* clone() const { return new BremsIntegral(*this); }

        // ----------------------------------------------------------------- //
        // Public methods
        // ----------------------------------------------------------------- //

        double CalculatedEdx(double energy);
};

} /* PROPOSAL */
