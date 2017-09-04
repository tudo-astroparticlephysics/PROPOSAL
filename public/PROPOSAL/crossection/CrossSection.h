
#pragma once

#include <vector>

#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Components.h"

namespace PROPOSAL
{

class Parametrization;

class CrossSection
{
    public:
        CrossSection(Parametrization&);
        CrossSection(const CrossSection&);
        virtual ~CrossSection();

        virtual CrossSection* clone() const = 0;

        // ----------------------------------------------------------------- //
        // Public methods
        // ----------------------------------------------------------------- //

        virtual double CalculatedEdx(double energy) = 0;
        virtual double CalculatedNdx(double energy) = 0;
        virtual double CalculatedNdx(double energy, double rnd) = 0;
        virtual double CalculateStochasticLoss(double energy, double rnd1, double rnd2) = 0;

    protected:
        typedef std::vector<Integral> IntegralVec;
        typedef std::vector<Components::Component*> ComponentVec;

        // ----------------------------------------------------------------- //
        // Protected methods
        // ----------------------------------------------------------------- //

        virtual double CalculateStochasticLoss(double energy, double rnd1) = 0;

        // ----------------------------------------------------------------- //
        // Protected member
        // ----------------------------------------------------------------- //

        Parametrization& parametrization;

        std::vector<double> prob_for_component_; //!< probability for each medium component to interact with the particle (formerly h_)
        double sum_of_rates_;

        Integral dedx_integral_;
        IntegralVec  dndx_integral_;

        const ComponentVec& components;

        double rnd_;    //!< This random number will be stored in CalculateDNdx to avoid calculate dNdx a second time in ClaculateSochasticLoss when it is already done
};

} /* PROPOSAL */
