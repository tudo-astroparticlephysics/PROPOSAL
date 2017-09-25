
#pragma once


#include "PROPOSAL/crossection/CrossSectionInterpolant.h"

namespace PROPOSAL
{

class Bremsstrahlung;

class BremsInterpolant: public CrossSectionInterpolant
{
    public:
        BremsInterpolant(const Bremsstrahlung&, InterpolationDef);
        BremsInterpolant(const BremsInterpolant&);
        virtual ~BremsInterpolant();

        CrossSection* clone() const { return new BremsInterpolant(*this); }

        // ----------------------------------------------------------------- //
        // Public methods
        // ----------------------------------------------------------------- //

        double CalculatedEdx(double energy);
};

} /* PROPOSAL */
