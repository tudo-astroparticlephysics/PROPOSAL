#pragma once


#include "PROPOSAL/crossection/CrossSectionInterpolant.h"

namespace PROPOSAL
{

class EpairInterpolant: public CrossSectionInterpolant
{
    public:
        EpairInterpolant(const Parametrization&, InterpolationDef = InterpolationDef());
        EpairInterpolant(const EpairInterpolant&);
        virtual ~EpairInterpolant();

        CrossSection* clone() const { return new EpairInterpolant(*this); }

        // ----------------------------------------------------------------- //
        // Public methods
        // ----------------------------------------------------------------- //

        double CalculatedEdx(double energy);
};

} /* PROPOSAL */
