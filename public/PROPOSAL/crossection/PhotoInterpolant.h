#pragma once


#include "PROPOSAL/crossection/CrossSectionInterpolant.h"

namespace PROPOSAL
{

class Photonuclear;

class PhotoInterpolant: public CrossSectionInterpolant
{
    public:
        PhotoInterpolant(const Photonuclear&, InterpolationDef);
        PhotoInterpolant(const PhotoInterpolant&);
        virtual ~PhotoInterpolant();

        CrossSection* clone() const { return new PhotoInterpolant(*this); }

        // ----------------------------------------------------------------- //
        // Public methods
        // ----------------------------------------------------------------- //

        double CalculatedEdx(double energy);
};

} /* PROPOSAL */
