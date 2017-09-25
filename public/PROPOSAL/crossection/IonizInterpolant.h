
#pragma once


#include "PROPOSAL/crossection/CrossSectionInterpolant.h"

namespace PROPOSAL
{

class Ionization;

class IonizInterpolant: public CrossSectionInterpolant
{
    public:
        IonizInterpolant(const Ionization&, InterpolationDef);
        IonizInterpolant(const IonizInterpolant&);
        virtual ~IonizInterpolant();

        CrossSection* clone() const { return new IonizInterpolant(*this); }

        // ----------------------------------------------------------------- //
        // Public methods
        // ----------------------------------------------------------------- //

        double CalculatedEdx(double energy);
        virtual double CalculatedNdx(double energy);
        virtual double CalculatedNdx(double energy, double rnd);

        // Needed to initialize interpolation
        double FunctionToBuildDNdxInterpolant(double energy, int component);
        virtual double FunctionToBuildDNdxInterpolant2D(double energy, double v, Integral&, int component);

    private:
        virtual double CalculateStochasticLoss(double energy, double rnd1);
        virtual void InitdNdxInerpolation(const InterpolationDef& def);
};

} /* PROPOSAL */
