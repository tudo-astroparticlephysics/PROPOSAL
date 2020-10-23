
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
#include <memory>

namespace PROPOSAL {
class Interpolant;
} // namespace PROPOSAL

namespace PROPOSAL {
struct InterpolantBuilder {
    InterpolantBuilder() = default;
    virtual ~InterpolantBuilder() = default;

    virtual std::unique_ptr<Interpolant> build() const = 0;
};

struct Interpolant1DBuilder : public InterpolantBuilder {
    struct Definition {
        int nodes = 100;
        int romberg = 5;
        int rombergY = 5;
        bool logSubst = false;
        bool isLog = true;
        bool rational = false;
        bool rationalY = false;
        bool relative = false;
        bool relativeY = false;
        double xmin = 0;
        double xmax = 1.e14;

        std::function<double(double)> function1d = nullptr;

        size_t GetHash() ;
    };

    Interpolant1DBuilder() = default;
    Interpolant1DBuilder(Interpolant1DBuilder::Definition);

    std::unique_ptr<Interpolant> build() const override;

private:
    std::unique_ptr<Definition> def;
};

struct Interpolant2DBuilder : public Interpolant1DBuilder {

    struct Definition : public Interpolant1DBuilder::Definition {
        int nodes2 = 100;
        int romberg2 = 5;
        bool isLog2 = false;
        bool rational2 = false;
        bool relative2 = false;
        double x2min = 0.0;
        double x2max = 1.0;

        std::function<double(double, double)> function2d = nullptr;

        size_t GetHash() ;
    };

    Interpolant2DBuilder() = default;
    Interpolant2DBuilder(Interpolant2DBuilder::Definition def);

    std::unique_ptr<Interpolant> build() const final;

private:
    std::unique_ptr<Definition> def;
};

} // namespace PROPOSAL
