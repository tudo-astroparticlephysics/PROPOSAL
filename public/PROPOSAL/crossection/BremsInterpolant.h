
#pragma once


#include "PROPOSAL/crossection/CrossSectionInterpolant.h"

namespace PROPOSAL
{

class BremsInterpolant: public CrossSectionInterpolant
{
    public:
        BremsInterpolant(Parametrization&);
        BremsInterpolant(const BremsInterpolant&);
        virtual ~BremsInterpolant();

        CrossSection* clone() const { return new BremsInterpolant(*this); }

        // ----------------------------------------------------------------- //
        // Public methods
        // ----------------------------------------------------------------- //

        double CalculatedEdx(double energy);
};

} /* PROPOSAL */
