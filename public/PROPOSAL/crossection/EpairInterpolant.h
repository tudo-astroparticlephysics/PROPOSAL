#pragma once


#include "PROPOSAL/crossection/CrossSectionInterpolant.h"

namespace PROPOSAL
{

class EpairProduction;

class EpairInterpolant: public CrossSectionInterpolant
{
    public:
        EpairInterpolant(const EpairProduction&, InterpolationDef);
        EpairInterpolant(const EpairInterpolant&);
        virtual ~EpairInterpolant();

        CrossSection* clone() const { return new EpairInterpolant(*this); }

        // ----------------------------------------------------------------- //
        // Public methods
        // ----------------------------------------------------------------- //

        double CalculatedEdx(double energy);
};

} /* PROPOSAL */
