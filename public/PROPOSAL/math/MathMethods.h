
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
#include <iostream>
#include "PROPOSAL/Logging.h"
#include <vector>

namespace PROPOSAL {

/// @brief Netwon-Raphson method and bisection to find the root of the function f
/// @param f    Function to find the root. Root must exist and be inside the interval [x1, x2]
/// @param df   Derivative of f
/// @param x1   lower limit of x to find the root
/// @param x2   upper limit of x to find the root
/// @param xacc convergence criterion: if $ \frac{f(x)}{df(x)} < xacc $ than accept x as the root
/// @return root of x

double NewtonRaphson(std::function<double(double)> f, std::function<double(double)> df, double x1, double x2, double xacc = 1.e-6);

struct SplineCoefficients{
    SplineCoefficients()
        : _a(0)
        , _b(0)
        , _c(0)
        , _d(0)
        , _x_t(0)
        {
        }

    SplineCoefficients(double a, double b, double c, double d, double x_t)
        : _a(a)
        , _b(b)
        , _c(c)
        , _d(d)
        , _x_t(x_t)
        {
        }

        double _a;
        double _b;
        double _c;
        double _d;
        double _x_t;

        // Spline has the form S(x) = a + b * (x - x_t) + c * (x - x_t)**2 + d * (x - x_t)**3

    };


/// @brief Calculates the coefficients for the Natural Cubic Splines corresponding to the coordinates (x,y)
/// @param x x-coordinates
/// @param y y-coordinates
/// @return Coefficients of the splines i in the subinterval [x_i, x_i+1]

std::vector<SplineCoefficients> CalculateSpline(std::vector<double> x, std::vector<double> y);


std::pair<std::vector<double>, std::vector<double>> ParseSplineCoordinates(const std::string&);

} // namespace PROPOSAL
