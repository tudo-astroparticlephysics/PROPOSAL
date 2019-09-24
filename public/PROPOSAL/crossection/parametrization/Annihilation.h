
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

#include <functional>
#include <cmath>
#include <fstream>

#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/methods.h"


namespace PROPOSAL {

    class Interpolant;

    class Annihilation : public Parametrization
    {
    public:
        Annihilation(const ParticleDef&, const Medium&, double multiplier);
        Annihilation(const Annihilation&);
        virtual ~Annihilation();

        virtual Parametrization* clone() const = 0;

        // ----------------------------------------------------------------- //
        // Public methods
        // ----------------------------------------------------------------- //

        virtual double DifferentialCrossSection(double energy, double v) = 0;

        virtual IntegralLimits GetIntegralLimits(double energy);

        virtual size_t GetHash() const;

    protected:
        bool compare(const Parametrization&) const;

    };


    class AnnihilationHeitler : public Annihilation
    {
    public:
        typedef std::vector<Interpolant*> InterpolantVec;

        AnnihilationHeitler(const ParticleDef&, const Medium&, double multiplier);
        AnnihilationHeitler(const AnnihilationHeitler&);
        virtual ~AnnihilationHeitler();

        virtual Parametrization* clone() const { return new AnnihilationHeitler(*this); }
        static Annihilation* create(const ParticleDef& particle_def,
                                       const Medium& medium,
                                       double multiplier)
        {
            return new AnnihilationHeitler(particle_def, medium, multiplier);
        }

        virtual double DifferentialCrossSection(double energy, double v);

        const std::string& GetName() const { return name_; }

    protected:
        static const std::string name_;
    };






} // namespace PROPOSAL
