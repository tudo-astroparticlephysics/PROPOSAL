
#pragma once

#include "PROPOSAL/crossection/parametrization/bremsstrahlung/Bremsstrahlung.h"

namespace PROPOSAL {
class BremsAndreevBezrukovBugaev : public Bremsstrahlung
{
    public:
    BremsAndreevBezrukovBugaev(const ParticleDef&, const Medium&, const EnergyCutSettings&, Definition = Definition());
    BremsAndreevBezrukovBugaev(const BremsAndreevBezrukovBugaev&);
    virtual ~BremsAndreevBezrukovBugaev();

    virtual Parametrization* clone() const { return new BremsAndreevBezrukovBugaev(*this); }

    // ----------------------------------------------------------------- //
    // Public methods
    // ----------------------------------------------------------------- //

    virtual double CalculateParametrization(double energy, double v);
};
} /* PROPOSAL */
