
#pragma once

#include "PROPOSAL/crossection/parametrization/bremsstrahlung/Bremsstrahlung.h"

namespace PROPOSAL {
class BremsPetrukhinShestakov : public Bremsstrahlung
{
    public:
    BremsPetrukhinShestakov(const ParticleDef&, const Medium&, const EnergyCutSettings&, Definition = Definition());
    BremsPetrukhinShestakov(const BremsPetrukhinShestakov&);
    virtual ~BremsPetrukhinShestakov();

    virtual Parametrization* clone() const { return new BremsPetrukhinShestakov(*this); }

    // ----------------------------------------------------------------- //
    // Public methods
    // ----------------------------------------------------------------- //

    virtual double CalculateParametrization(double energy, double v);
};
} /* PROPOSAL */
