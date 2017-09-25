#pragma once


#include "PROPOSAL/crossection/CrossSectionIntegral.h"

namespace PROPOSAL
{

class EpairProduction;

class EpairIntegral: public CrossSectionIntegral
{
    public:
        EpairIntegral(const EpairProduction&);
        EpairIntegral(const EpairIntegral&);
        virtual ~EpairIntegral();

        CrossSection* clone() const { return new EpairIntegral(*this); }

        // ----------------------------------------------------------------- //
        // Public methods
        // ----------------------------------------------------------------- //

        double CalculatedEdx(double energy);

    private:
        // ----------------------------------------------------------------------------
        /// @brief calls FunctionToDEdxIntegral from the parametrization
        // ----------------------------------------------------------------------------
        double FunctionToDEdxIntegralReverse(double energy, double v);
};

} /* PROPOSAL */
