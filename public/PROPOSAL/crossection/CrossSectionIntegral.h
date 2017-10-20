
#pragma once

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/math/Integral.h"

namespace PROPOSAL
{

class CrossSectionIntegral: public CrossSection
{
    public:
        typedef std::vector<Integral> IntegralVec;

    public:
        CrossSectionIntegral(const DynamicData::Type&, const Parametrization&);
        CrossSectionIntegral(const CrossSectionIntegral&);
        virtual ~CrossSectionIntegral();

        virtual CrossSection* clone() const = 0;

        virtual double CalculatedEdx(double energy) = 0;
        virtual double CalculatedE2dx(double energy);
        virtual double CalculatedNdx(double energy);
        virtual double CalculatedNdx(double energy, double rnd);
        double CalculateStochasticLoss(double energy, double rnd1, double rnd2);

    protected:
        virtual bool compare(const CrossSection&) const;

        Integral dedx_integral_;
        Integral de2dx_integral_;
        IntegralVec  dndx_integral_;

        virtual double CalculateStochasticLoss(double energy, double rnd1);
};

} /* PROPOSAL */
