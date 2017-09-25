
#pragma once


#include "PROPOSAL/crossection/CrossSectionIntegral.h"

namespace PROPOSAL
{

class Ionization;

class IonizIntegral: public CrossSectionIntegral
{
    public:
        IonizIntegral(const Ionization&);
        IonizIntegral(const IonizIntegral&);
        virtual ~IonizIntegral();

        CrossSection* clone() const { return new IonizIntegral(*this); }

        // ----------------------------------------------------------------- //
        // Public methods
        // ----------------------------------------------------------------- //

        double CalculatedEdx(double energy);
        double CalculatedE2dx(double energy);
        double CalculatedNdx(double energy);
        double CalculatedNdx(double energy, double rnd);

    private:
        virtual double CalculateStochasticLoss(double energy, double rnd1);
        double Delta(double beta, double gamma);
};

} /* PROPOSAL */
