
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

#include "PROPOSAL/crossection/CrossSectionInterpolant.h"

namespace PROPOSAL {

    class Compton;

    class ComptonInterpolant : public CrossSectionInterpolant
    {
    public:
        ComptonInterpolant(const Compton&, InterpolationDef);
        ComptonInterpolant(const ComptonInterpolant&);
        virtual ~ComptonInterpolant();

        CrossSection* clone() const { return new ComptonInterpolant(*this); }

        // ----------------------------------------------------------------- //
        // Public methods
        // ----------------------------------------------------------------- //

        double CalculatedEdx(double energy);
        virtual double FunctionToBuildDNdxInterpolant2D(double energy, double v, Integral&, int component);
        virtual double CalculateCumulativeCrossSection(double energy, int component, double v);
        virtual std::pair<double, double> StochasticDeflection(double energy, double energy_loss);

    private:
        virtual double CalculateStochasticLoss(double energy, double rnd1);
        virtual void InitdNdxInterpolation(const InterpolationDef& def);

    };

} // namespace PROPOSAL
