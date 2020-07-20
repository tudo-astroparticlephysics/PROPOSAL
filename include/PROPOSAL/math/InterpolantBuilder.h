
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
#include <type_traits>

using std::function;
using std::remove_reference;
using std::unique_ptr;

namespace PROPOSAL {
class Interpolant;
} // namespace PROPOSAL

namespace PROPOSAL {
class InterpolantBuilder {
public:
    InterpolantBuilder() = default;
    virtual ~InterpolantBuilder() = default;

    virtual unique_ptr<Interpolant> build() = 0;
};

class Interpolant1DBuilder : public InterpolantBuilder {
public:
    template <typename T> Interpolant1DBuilder(T&& def);

    struct Definition {
        function<double(double)> function1d = nullptr;
        int max = 200;
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
    };
    unique_ptr<Interpolant> build() override;

private:
    unique_ptr<Definition> builder_def;
};

class Interpolant2DBuilder : public InterpolantBuilder {

public:
    template <typename T> Interpolant2DBuilder(T&& def);

    struct Definition {
        function<double(double, double)> function2d = nullptr;
        int max1 = 1.0;
        double x1min = 1.0;
        double x1max = 1.0;
        int max2 = 1.0;
        double x2min = 1.0;
        double x2max = 1.0;
        int romberg1 = 1;
        bool rational1 = false;
        bool relative1 = false;
        bool isLog1 = false;
        int romberg2 = 1;
        bool rational2 = false;
        bool relative2 = false;
        bool isLog2 = false;
        int rombergY = 1;
        bool rationalY = false;
        bool relativeY = false;
        bool logSubst = false;
    };

    unique_ptr<Interpolant> build() override;

private:
    unique_ptr<Definition> builder_def;
};

template <typename T>
Interpolant1DBuilder::Interpolant1DBuilder(T&& def)
    : builder_def(new Definition(def))
{
}

template <typename T>
Interpolant2DBuilder::Interpolant2DBuilder(T&& def)
    : builder_def(new Definition(def))
{
}
} // namespace PROPOSAL
