
#pragma once

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/methods.h"

namespace PROPOSAL {

class Integral;
class Interpolant;

class CrossSectionInterpolant : public CrossSection
{
public:
    CrossSectionInterpolant(const DynamicData::Type&, const Parametrization&);
    CrossSectionInterpolant(const CrossSectionInterpolant&);
    virtual ~CrossSectionInterpolant();

    virtual CrossSection* clone() const = 0;

    virtual double CalculatedEdx(double energy) = 0;
    virtual double CalculatedE2dx(double energy);
    virtual double CalculatedNdx(double energy);
    virtual double CalculatedNdx(double energy, double rnd);
    virtual double CalculateStochasticLoss(double energy, double rnd1, double rnd2);

    // Needed to initialize interpolation
    virtual double FunctionToBuildDNdxInterpolant(double energy, int component);
    virtual double FunctionToBuildDNdxInterpolant2D(double energy, double v, Integral&, int component);

protected:
    virtual bool compare(const CrossSection&) const;

    typedef std::vector<Interpolant*> InterpolantVec;

    virtual double CalculateStochasticLoss(double energy, double rnd1);
    virtual void InitdNdxInerpolation(const InterpolationDef& def);

    Interpolant* dedx_interpolant_;
    Interpolant* de2dx_interpolant_;
    InterpolantVec dndx_interpolant_1d_; // Stochastic dNdx()
    InterpolantVec dndx_interpolant_2d_; // Stochastic dNdx()
};

} // namespace PROPOSAL
