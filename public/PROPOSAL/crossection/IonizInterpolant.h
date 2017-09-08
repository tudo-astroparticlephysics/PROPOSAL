
#pragma once


#include "PROPOSAL/crossection/CrossSectionInterpolant.h"

namespace PROPOSAL
{

class IonizInterpolant: public CrossSectionInterpolant
{
    public:
        IonizInterpolant(Parametrization&);
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
        // double FunctionToBuildDEdxInterpolant(double energy);
        double FunctionToBuildDNdxInterpolant(double energy, int component);
        double FunctionToBuildDNdxInterpolant2D(double energy, double v, int component);

    private:
        virtual double CalculateStochasticLoss(double energy, double rnd1);
        //TODO(mario): Remove here Thu 2017/09/07
        // double Delta(double beta, double gamma);
};

} /* PROPOSAL */
