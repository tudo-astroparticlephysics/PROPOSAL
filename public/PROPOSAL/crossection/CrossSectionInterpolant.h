
#pragma once

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/math/Interpolant.h"

namespace PROPOSAL
{

class CrossSectionInterpolant: public CrossSection
{
    public:
        CrossSectionInterpolant(Parametrization&);
        CrossSectionInterpolant(const CrossSectionInterpolant&);
        virtual ~CrossSectionInterpolant();

        virtual CrossSection* clone() const = 0;

        virtual double CalculatedEdx(double energy) = 0;
        virtual double CalculatedNdx(double energy);
        virtual double CalculatedNdx(double energy, double rnd);
        virtual double CalculateStochasticLoss(double energy, double rnd1, double rnd2);

        // Needed to initialize interpolation
        // virtual double FunctionToBuildDEdxInterpolant(double energy) = 0;
        virtual double FunctionToBuildDNdxInterpolant(double energy, int component);
        virtual double FunctionToBuildDNdxInterpolant2D(double energy, double v, int component);

    protected:
        typedef std::vector<Interpolant*> InterpolantVec;

        virtual double CalculateStochasticLoss(double energy, double rnd1);
        // void InitInterpolation(std::string filename, bool raw);

        Interpolant* dedx_interpolant_;
        InterpolantVec dndx_interpolant_1d_; //Stochastic dNdx()
        InterpolantVec dndx_interpolant_2d_; //Stochastic dNdx()
};

} /* PROPOSAL */
