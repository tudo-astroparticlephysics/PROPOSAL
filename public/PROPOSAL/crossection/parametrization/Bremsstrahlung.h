
/******************************************************************************
 *                                                                            *
 * This file is part of the simulation tool PROPOSAL.                         *
 *                                                                            *
 * Copyright (C) 2017 TU Dortmund University, Department of Physics,          *
 *                    Chair Experimental Physics 5b                           *
 *                                                                            *
 * This software may be modified and distributed under the terms of a         *
 * modified GNU Lesser General Public Licence version 3 (LGPL),               *
 * copied verbatim in the file "LICENSE".                                     *
 *                                                                            *
 * Modifcations to the LGPL License:                                          *
 *                                                                            *
 *      1. The user shall acknowledge the use of PROPOSAL by citing the       *
 *         following reference:                                               *
 *                                                                            *
 *         J.H. Koehne et al.  Comput.Phys.Commun. 184 (2013) 2070-2090 DOI:  *
 *         10.1016/j.cpc.2013.04.001                                          *
 *                                                                            *
 *      2. The user should report any bugs/errors or improvments to the       *
 *         current maintainer of PROPOSAL or open an issue on the             *
 *         GitHub webpage                                                     *
 *                                                                            *
 *         "https://github.com/tudo-astroparticlephysics/PROPOSAL"            *
 *                                                                            *
 ******************************************************************************/


#pragma once

#include "PROPOSAL/crossection/parametrization/Parametrization.h"

#define BREMSSTRAHLUNG_DEF(param)                                                                                      \
    class Brems##param : public Bremsstrahlung                                                                         \
    {                                                                                                                  \
    public:                                                                                                            \
        Brems##param(const ParticleDef&, const Medium&, const EnergyCutSettings&, double multiplier, bool lpm);        \
        Brems##param(const Brems##param&);                                                                             \
        ~Brems##param();                                                                                               \
                                                                                                                       \
        Parametrization* clone() const { return new Brems##param(*this); }                                             \
        static Bremsstrahlung* create(const ParticleDef& particle_def,                                                 \
                                      const Medium& medium,                                                            \
                                      const EnergyCutSettings& cuts,                                                   \
                                      double multiplier,                                                               \
                                      bool lpm)                                                                        \
        {                                                                                                              \
            return new Brems##param(particle_def, medium, cuts, multiplier, lpm);                                      \
        }                                                                                                              \
                                                                                                                       \
        double CalculateParametrization(double energy, double v);                                                      \
                                                                                                                       \
        const std::string& GetName() const { return name_; }                                                           \
                                                                                                                       \
    private:                                                                                                           \
        static const std::string name_;                                                                                \
    };

namespace PROPOSAL {
class Bremsstrahlung : public Parametrization
{
public:
    Bremsstrahlung(const ParticleDef&, const Medium&, const EnergyCutSettings&, double multiplier, bool lpm);
    Bremsstrahlung(const Bremsstrahlung&);
    virtual ~Bremsstrahlung();

    virtual Parametrization* clone() const = 0;

    // ----------------------------------------------------------------- //
    // Public methods
    // ----------------------------------------------------------------- //

    virtual double DifferentialCrossSection(double energy, double v);
    virtual double CalculateParametrization(double energy, double v) = 0;

    virtual IntegralLimits GetIntegralLimits(double energy);

    // ----------------------------------------------------------------- //
    // Getter
    // ----------------------------------------------------------------- //

    virtual size_t GetHash() const;

protected:
    virtual bool compare(const Parametrization&) const;
    virtual void print(std::ostream&) const;

    // ----------------------------------------------------------------- //
    // Protected methods
    // ----------------------------------------------------------------- //

    double lpm(double energy, double v);

    // ----------------------------------------------------------------- //
    // Protected member
    // ----------------------------------------------------------------- //

    bool lorenz_;       /// enable lorenz cut
    double lorenz_cut_; /// in [MeV] // - set to 1.e6 in Constructor
    bool init_lpm_effect_;
    bool lpm_;
    double eLpm_;
};

// ------------------------------------------------------------------------- //
// Declare the specific parametrizations
// ------------------------------------------------------------------------- //

BREMSSTRAHLUNG_DEF(PetrukhinShestakov)
BREMSSTRAHLUNG_DEF(KelnerKokoulinPetrukhin)
BREMSSTRAHLUNG_DEF(CompleteScreening)
BREMSSTRAHLUNG_DEF(AndreevBezrukovBugaev)
BREMSSTRAHLUNG_DEF(SandrockSoedingreksoRhode)

#undef BREMSSTRAHLUNG_DEF

} // namespace PROPOSAL
