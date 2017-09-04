
#pragma once

#include "PROPOSAL/crossection/parametrization/bremsstrahlung/Bremsstrahlung.h"

namespace PROPOSAL {
class BremsKelnerKokoulinPetrukhin : public Bremsstrahlung
{
    public:
    BremsKelnerKokoulinPetrukhin(const ParticleDef&, const Medium&, const EnergyCutSettings&, Definition = Definition());
    BremsKelnerKokoulinPetrukhin(const BremsKelnerKokoulinPetrukhin&);
    virtual ~BremsKelnerKokoulinPetrukhin();

    virtual Parametrization* clone() const { return new BremsKelnerKokoulinPetrukhin(*this); }

    // ----------------------------------------------------------------- //
    // Public methods
    // ----------------------------------------------------------------- //

    virtual double CalculateParametrization(double energy, double v);
};
} /* PROPOSAL */
