
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

namespace PROPOSAL {
class Ionization : public Parametrization
{
public:
    Ionization(const ParticleDef&, std::shared_ptr<const Medium>, const EnergyCutSettings&, double multiplier);
    Ionization(const Ionization&);
    virtual ~Ionization();

    virtual Parametrization* clone() const = 0;

    // ----------------------------------------------------------------- //
    // Public methods
    // ----------------------------------------------------------------- //
    double Delta(double beta, double gamma);
    double DifferentialCrossSection(double energy, double v) = 0;
private:

};

// ----------------------------------------------------------------- //
// Spefific Parametrization
// ----------------------------------------------------------------- //

class IonizBetheBlochRossi : public Ionization
{
public:
    IonizBetheBlochRossi(const ParticleDef&, std::shared_ptr<const Medium>, const EnergyCutSettings&, double multiplier);
    IonizBetheBlochRossi(const IonizBetheBlochRossi&);
    ~IonizBetheBlochRossi();

    Parametrization* clone() const { return new IonizBetheBlochRossi(*this); }
    static Ionization* create(const ParticleDef& particle_def,
                                      std::shared_ptr<const Medium> medium,
                                      const EnergyCutSettings& cuts,
                                      double multiplier)
        {
            return new IonizBetheBlochRossi(particle_def, medium, cuts, multiplier);
        }

     IntegralLimits GetIntegralLimits(double energy);
     double DifferentialCrossSection(double energy, double v);
     double FunctionToDEdxIntegral(double energy, double v);

     const std::string& GetName() const { return name_; }

private:
     double InelCorrection(double energy, double v);
     double CrossSectionWithoutInelasticCorrection(double energy, double v);
     static const std::string name_;
};

class IonizBergerSeltzerBhabha : public Ionization
{
public:
    IonizBergerSeltzerBhabha(const ParticleDef&, std::shared_ptr<const Medium>, const EnergyCutSettings&, double multiplier);
    IonizBergerSeltzerBhabha(const IonizBergerSeltzerBhabha&);
    ~IonizBergerSeltzerBhabha();

    Parametrization* clone() const { return new IonizBergerSeltzerBhabha(*this); }
    static Ionization* create(const ParticleDef& particle_def,
                                    std::shared_ptr<const Medium> medium,
                                    const EnergyCutSettings& cuts,
                                    double multiplier)
        {
            return new IonizBergerSeltzerBhabha(particle_def, medium, cuts, multiplier);
        }

    IntegralLimits GetIntegralLimits(double energy);
    double DifferentialCrossSection(double energy, double v);
    double FunctionToDEdxIntegral(double energy, double v);

    const std::string& GetName() const { return name_; }

private:
    static const std::string name_;
};

class IonizBergerSeltzerMoller : public Ionization
    {
    public:
    IonizBergerSeltzerMoller(const ParticleDef&, std::shared_ptr<const Medium>, const EnergyCutSettings&, double multiplier);
    IonizBergerSeltzerMoller(const IonizBergerSeltzerMoller&);
    ~IonizBergerSeltzerMoller();

    Parametrization* clone() const { return new IonizBergerSeltzerMoller(*this); }
        static Ionization* create(const ParticleDef& particle_def,
                                  std::shared_ptr<const Medium> medium,
                                  const EnergyCutSettings& cuts,
                                  double multiplier)
        {
            return new IonizBergerSeltzerMoller(particle_def, medium, cuts, multiplier);
        }

        IntegralLimits GetIntegralLimits(double energy);
        double DifferentialCrossSection(double energy, double v);
        double FunctionToDEdxIntegral(double energy, double v);

        const std::string& GetName() const { return name_; }

    private:
        static const std::string name_;
    };

} // namespace PROPOSAL
