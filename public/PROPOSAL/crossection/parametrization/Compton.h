
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
    class Compton : public Parametrization
    {
    public:
        Compton(const ParticleDef&, std::shared_ptr<const Medium>, const EnergyCutSettings&, double multiplier);
        Compton(const Compton&);
        virtual ~Compton();

        virtual Parametrization* clone() const = 0;

        // ----------------------------------------------------------------- //
        // Public methods
        // ----------------------------------------------------------------- //

        virtual double DifferentialCrossSection(double energy, double v) = 0;

        virtual IntegralLimits GetIntegralLimits(double energy);

        // ----------------------------------------------------------------- //
        // Getter
        // ----------------------------------------------------------------- //

        virtual size_t GetHash() const;

    protected:
        virtual bool compare(const Parametrization&) const;

    };

// ------------------------------------------------------------------------- //
// Declare the specific parametrizations
// ------------------------------------------------------------------------- //


    class ComptonKleinNishina : public Compton
    {
    public:
        ComptonKleinNishina(const ParticleDef&, std::shared_ptr<const Medium>, const EnergyCutSettings&, double multiplier);
        ComptonKleinNishina(const ComptonKleinNishina&);
        ~ComptonKleinNishina();

        Parametrization* clone() const { return new ComptonKleinNishina(*this); }
        static Compton* create(const ParticleDef& particle_def,
                                      std::shared_ptr<const Medium> medium,
                                      const EnergyCutSettings& cuts,
                                      double multiplier)
        {
            return new ComptonKleinNishina(particle_def, medium, cuts, multiplier);
        }

        double DifferentialCrossSection(double energy, double v);

        const std::string& GetName() const { return name_; }

    private:
        virtual bool compare(const Parametrization&) const;

        static const std::string name_;
    };


} // namespace PROPOSAL
