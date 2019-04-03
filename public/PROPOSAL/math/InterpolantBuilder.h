
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
#include <vector>

namespace PROPOSAL {

class Interpolant;

// ----------------------------------------------------------------------------
/// @brief Base class for interpolantBuilder
///
/// The Builder is used for initializing the interpolants.
// ----------------------------------------------------------------------------
class InterpolantBuilder
{
public:
    // Default values for variables
    static const int default_max;
    static const double default_xmin;
    static const double default_xmax;

    static const int default_romberg;
    static const bool default_rational;
    static const bool default_relative;
    static const bool default_isLog;

    static const int default_rombergY;
    static const bool default_rationalY;
    static const bool default_relativeY;
    static const bool default_logSubst;

    InterpolantBuilder() {}
    virtual ~InterpolantBuilder() {}

    virtual Interpolant* build() = 0;
};

// ----------------------------------------------------------------------------
/// @brief
// ----------------------------------------------------------------------------
class Interpolant1DBuilder : public InterpolantBuilder
{
public:
    typedef std::function<double(double)> Function1D;
    static const Function1D default_function1d;

public:
    // Constructor
    Interpolant1DBuilder();
    Interpolant1DBuilder(const Interpolant1DBuilder&);

    // Setter to build the Interpolant
    Interpolant1DBuilder& SetFunction1D(Function1D val)
    {
        function1d = val;
        return *this;
    }

    Interpolant1DBuilder& SetMax(const int val)
    {
        max = val;
        return *this;
    }
    Interpolant1DBuilder& SetXMin(const double val)
    {
        xmin = val;
        return *this;
    }
    Interpolant1DBuilder& SetXMax(const double val)
    {
        xmax = val;
        return *this;
    }

    Interpolant1DBuilder& SetRomberg(const int val)
    {
        romberg = val;
        return *this;
    }
    Interpolant1DBuilder& SetRational(const bool val)
    {
        rational = val;
        return *this;
    }
    Interpolant1DBuilder& SetRelative(const bool val)
    {
        relative = val;
        return *this;
    }
    Interpolant1DBuilder& SetIsLog(const bool val)
    {
        isLog = val;
        return *this;
    }

    Interpolant1DBuilder& SetRombergY(const int val)
    {
        rombergY = val;
        return *this;
    }
    Interpolant1DBuilder& SetRationalY(const bool val)
    {
        rationalY = val;
        return *this;
    }
    Interpolant1DBuilder& SetRelativeY(const bool val)
    {
        relativeY = val;
        return *this;
    }
    Interpolant1DBuilder& SetLogSubst(const bool val)
    {
        logSubst = val;
        return *this;
    }

    // prepare specific frequently desired Product
    // returns Builder for shorthand inline usage (same way as cout <<)
    // Builder& setProductP(){
    // 	this->i = 42;
    // 	this->f = -1.0f/12.0f;
    // 	this->c = '@';
    //
    // 	return *this;
    // }

    Interpolant* build();

private:
    Function1D function1d;

    int max;
    double xmin, xmax;

    int romberg;
    bool rational, relative, isLog;

    int rombergY;
    bool rationalY, relativeY, logSubst;
};

// ----------------------------------------------------------------------------
/// @brief
// ----------------------------------------------------------------------------
class Interpolant2DBuilder : public InterpolantBuilder
{
public:
    typedef std::function<double(double, double)> Function2D;
    static const Function2D default_function2d;

    // Constructor
    Interpolant2DBuilder();
    Interpolant2DBuilder(const Interpolant2DBuilder&);

    // Setter to build the Interpolant
    Interpolant2DBuilder& SetFunction2D(Function2D val)
    {
        function2d = val;
        return *this;
    }

    Interpolant2DBuilder& SetMax1(const int val)
    {
        max1 = val;
        return *this;
    }
    Interpolant2DBuilder& SetX1Min(const double val)
    {
        x1min = val;
        return *this;
    }
    Interpolant2DBuilder& SetX1Max(const double val)
    {
        x1max = val;
        return *this;
    }

    Interpolant2DBuilder& SetMax2(const int val)
    {
        max2 = val;
        return *this;
    }
    Interpolant2DBuilder& SetX2Min(const double val)
    {
        x2min = val;
        return *this;
    }
    Interpolant2DBuilder& SetX2Max(const double val)
    {
        x2max = val;
        return *this;
    }

    Interpolant2DBuilder& SetRomberg1(const int val)
    {
        romberg1 = val;
        return *this;
    }
    Interpolant2DBuilder& SetRational1(const bool val)
    {
        rational1 = val;
        return *this;
    }
    Interpolant2DBuilder& SetRelative1(const bool val)
    {
        relative1 = val;
        return *this;
    }

    Interpolant2DBuilder& SetIsLog1(const bool val)
    {
        isLog1 = val;
        return *this;
    }

    Interpolant2DBuilder& SetRomberg2(const int val)
    {
        romberg2 = val;
        return *this;
    }
    Interpolant2DBuilder& SetRational2(const bool val)
    {
        rational2 = val;
        return *this;
    }
    Interpolant2DBuilder& SetRelative2(const bool val)
    {
        relative2 = val;
        return *this;
    }

    Interpolant2DBuilder& SetIsLog2(const bool val)
    {
        isLog2 = val;
        return *this;
    }

    Interpolant2DBuilder& SetRombergY(const int val)
    {
        rombergY = val;
        return *this;
    }
    Interpolant2DBuilder& SetRationalY(const bool val)
    {
        rationalY = val;
        return *this;
    }
    Interpolant2DBuilder& SetRelativeY(const bool val)
    {
        relativeY = val;
        return *this;
    }
    Interpolant2DBuilder& SetLogSubst(const bool val)
    {
        logSubst = val;
        return *this;
    }

    // prepare specific frequently desired Product
    // returns Builder for shorthand inline usage (same way as cout <<)
    // Builder& setProductP(){
    // 	this->i = 42;
    // 	this->f = -1.0f/12.0f;
    // 	this->c = '@';
    //
    // 	return *this;
    // }

    Interpolant* build();

private:
    Function2D function2d;

    int max1;
    double x1min, x1max;

    int max2;
    double x2min, x2max;

    int romberg1;
    bool rational1, relative1, isLog1;

    int romberg2;
    bool rational2, relative2, isLog2;

    int rombergY;
    bool rationalY, relativeY, logSubst;
};

// ----------------------------------------------------------------------------
/// @brief
// ----------------------------------------------------------------------------
    class Interpolant2DBuilder_array_as : public InterpolantBuilder
    {
    public:


        // Constructor
        Interpolant2DBuilder_array_as();
        Interpolant2DBuilder_array_as(const Interpolant2DBuilder_array_as&);

        // Setter to build the Interpolant
        Interpolant2DBuilder_array_as& Setx1(std::vector<double> vec)
        {
            x1 = vec;
            return *this;
        }

        Interpolant2DBuilder_array_as& Setx2(std::vector< std::vector<double> > vec)
        {
            x2 = vec;
            return *this;
        }

        Interpolant2DBuilder_array_as& Sety(std::vector< std::vector<double> > vec)
        {
            y = vec;
            return *this;
        }

        Interpolant2DBuilder_array_as& SetRomberg1(const int val)
        {
            romberg1 = val;
            return *this;
        }
        Interpolant2DBuilder_array_as& SetRational1(const bool val)
        {
            rational1 = val;
            return *this;
        }
        Interpolant2DBuilder_array_as& SetRelative1(const bool val)
        {
            relative1 = val;
            return *this;
        }

        Interpolant2DBuilder_array_as& SetRomberg2(const int val)
        {
            romberg2 = val;
            return *this;
        }
        Interpolant2DBuilder_array_as& SetRational2(const bool val)
        {
            rational2 = val;
            return *this;
        }
        Interpolant2DBuilder_array_as& SetRelative2(const bool val)
        {
            relative2 = val;
            return *this;
        }


        Interpolant* build();

    private:

        std::vector<double> x1;
        std::vector< std::vector<double> > x2;
        std::vector< std::vector<double> > y;

        int romberg1;
        bool rational1, relative1;

        int romberg2;
        bool rational2, relative2;

    };

} // namespace PROPOSAL
