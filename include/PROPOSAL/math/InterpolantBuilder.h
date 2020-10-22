
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
class InterpolantBuilder {
public:
    InterpolantBuilder() = default;
    virtual ~InterpolantBuilder() = default;

    virtual std::unique_ptr<Interpolant> build() const = 0;
};

struct Interpolant1DBuilder : public InterpolantBuilder {
    struct Definition {
        std::function<double(double)> function1d = nullptr;
        int nodes = 200;
        double xmin = 0;
        double xmax = 1.e14;
        int romberg = 5;
        bool rational = false;
        bool relative = false;
        bool isLog = true;
        int rombergY = 5;
        bool rationalY = false;
        bool relativeY = false;
        bool logSubst = false;

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
        std::function<double(double, double)> function2d = nullptr;
        int nodes2 = 1.0;
        double x2min = 1.0;
        double x2max = 1.0;
        int romberg2 = 1;
        bool rational2 = false;
        bool relative2 = false;
        bool isLog2 = false;

        size_t GetHash() ;
    };

    Interpolant2DBuilder() = default;
    Interpolant2DBuilder(Interpolant2DBuilder::Definition def);

    std::unique_ptr<Interpolant> build() const final;

private:
    std::unique_ptr<Definition> def;
};

} // namespace PROPOSAL
