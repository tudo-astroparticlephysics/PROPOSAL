
#pragma once

#include "PROPOSAL/crossection/parametrization/bremsstrahlung/Bremsstrahlung.h"

namespace PROPOSAL {
class BremsCompleteScreening : public Bremsstrahlung
{
    public:
    BremsCompleteScreening(const ParticleDef&, const Medium&, const EnergyCutSettings&, Definition = Definition());
    BremsCompleteScreening(const BremsCompleteScreening&);
    virtual ~BremsCompleteScreening();

    virtual Parametrization* clone() const { return new BremsCompleteScreening(*this); }

    // ----------------------------------------------------------------- //
    // Public methods
    // ----------------------------------------------------------------- //

    virtual double CalculateParametrization(double energy, double v);
};
} /* PROPOSAL */
