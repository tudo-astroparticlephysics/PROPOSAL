
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

#include <vector>
#include <string>
#include <fstream>

#include <functional>

namespace PROPOSAL {

/**
 *\class Interpolant
 *
 * This class provides routines for function interpolation. Include the function to be interpolated
 * in a class that implements the interface FunctionInt or FunctionInt2 (defined below).
 * Methods contained here are based on the Numerical Recipes (W. H. Press et al.).
 * <pre>
 * interface FunctionInt{
 *     double functionInt(double x);
 * }
 *
 * interface FunctionInt2{
 *     double functionInt(double x1, double x2);
 * }
 * </pre>
 * @author Dmitry Chirkin
 */

class Interpolant
{ /// implements FunctionInt{

private:
    const static double bigNumber_;
    const static double aBigNumber_;

    int romberg_, rombergY_;

    std::vector<double> iX_;
    std::vector<double> iY_;

    std::vector<std::vector<double> > iY2_;

    std::vector<double> c_;
    std::vector<double> d_;

    int max_;
    double xmin_, xmax_, step_;
    bool rational_, relative_;

    std::function<double(double)> function1d_;
    std::function<double(double, double)> function2d_;
    std::vector<Interpolant*> Interpolant_;

    int row_, starti_;
    bool rationalY_, relativeY_;

    bool reverse_, self_, flag_; // Self is setted to true in constructor
    bool isLog_, logSubst_;

    double precision_, worstX_;
    double precision2_, worstX2_;
    double precisionY_, worstY_;

    bool fast_; // Is setted to true in constructor

    double x_save_, y_save_; // Is setted to 1 and 0 in constructor

    //----------------------------------------------------------------------------//
    // Memberfunctions

    /*!
     * interpolates f(x) based on the values iY[i]=f(iX[i]) in the romberg-vicinity of x
     *
     * \param   x        position of the function
     * \param   start    start position of the sampling points for interpolation
     * \return  Interpolation result
     */
    double Interpolate(double x, int start);

    //----------------------------------------------------------------------------//

    /**
     * Exp(x) with CutOff.
     *
     * if x > exp(aBigNumber): exp(x) \n
     * else: 0
     *
     * \param    x
     * \return   exp(x) OR 0;
     */

    double Exp(double x);

    //----------------------------------------------------------------------------//

    /**
     * Log if not zero.
     *
     * if x > 0: log(x) \n
     * else: bigNumber
     *
     * \param    x
     * \return   log(x) OR bigNumber;
     */

    double Log(double x);

    //----------------------------------------------------------------------------//

    /**
     * Log it again.
     *
     * \param    x
     * \return   log(x) OR y_save;
     */

    double Slog(double x);

    //----------------------------------------------------------------------------//

    double Get2dFunctionFixedY(double x); // TOMSASZ

    //----------------------------------------------------------------------------//

    /**
     * Auxiliary class initializer.
     *
     * Initializes class for the 1-dimensional function.
     *
     * \param   max          number of sampling points
     * \param   xmin         lower limit of the interpolation table
     * \param   xmax         upper limit of the interpolation table
     * \param   romberg      order of interpolation
     * \param   rational     interpolate with rational function
     * \param   relative     save error relative to the function value
     * \param   isLog        substitute x = log(x)
     * \param   rombergY     order of interpolation for inverse interpolation(Find Limit)
     * \param   rationalY    interpolate with rational function
     * \param   relativeY    save error relative to the interpolated x-value
     * \param   logSubst     substitute f(x) = log(f(x))
     * \return
     */
    void InitInterpolant(int max,
                         double xmin,
                         double xmax,
                         int romberg,
                         bool rational,
                         bool relative,
                         bool isLog,
                         int rombergY,
                         bool rationalY,
                         bool relativeY,
                         bool logSubst);

    //----------------------------------------------------------------------------//

    /**
     * Defines a function for every row.
     *
     * \param    x
     * \return   Function value f(x)
     */

    double FunctionInt(double x);

    //----------------------------------------------------------------------------//

public:
    Interpolant(const Interpolant&);
    Interpolant& operator=(const Interpolant&);
    bool operator==(const Interpolant& interpolant) const;
    bool operator!=(const Interpolant& interpolant) const;
    /*!
     * Default Constructor.
     * Default Constructor, which does nothing.
     */
    Interpolant();

    /*!
     * Main constructor for the 1-dimensional functions.
     *
     * Main constructor and initializes class for the 1-dimensional function.
     *
     * \param   max          number of sampling points
     * \param   xmin         lower limit of the interpolation table
     * \param   xmax         upper limit of the interpolation table
     * \param   function2int function which will be interpolated
     * \param   romberg      Order of interpolation
     * \param   rational     interpolate with rational function?
     * \param   relative     save error relative to the function value?
     * \param   isLog        substitute x = log(x)?
     * \param   rombergY     Order of interpolation for inverse interpolation(Find Limit)
     * \param   rationalY    interpolate with rational function?
     * \param   relativeY    save error relative to the interpolated x-value?
     * \param   logSubst     substitute f(x) = log(f(x))?
     * \return
     */
    Interpolant(int max,
                double xmin,
                double xmax,
                std::function<double(double)> function1d,
                int romberg,
                bool rational,
                bool relative,
                bool isLog,
                int rombergY,
                bool rationalY,
                bool relativeY,
                bool logSubst);

    //----------------------------------------------------------------------------//

    /*!
     * Main constructor for the 2-dimensional functions.
     *
     * Main constructor and initializes class for the 2-dimensional function.
     *
     * \param   max          number of sampling points
     * \param   xmin         lower limit of the interpolation table
     * \param   xmax         upper limit of the interpolation table
     * \param   function2int function which will be interpolated
     * \param   romberg1     Order of interpolation f(x,a)
     * \param   rational1    interpolate with rational function?
     * \param   relative1    save error relative to the function value?
     * \param   isLog1       substitute x1 = log(x2)
     * \param   romberg2     Order of interpolation for f(a,x)
     * \param   rational2    interpolate with rational function?
     * \param   relative2    save error relative to the function value?
     * \param   isLog2       substitute x2 = log(x2)
     * \param   rombergY     Order of interpolation for inverse interpolation(Find Limit)
     * \param   rationalY    interpolate with rational function?
     * \param   relativeY    save error relative to the interpolated x-value
     * \param   logSubst     substitute f(x1,x2) = log(f(x1,x2))
     * \return
     */
    Interpolant(int max1,
                double x1min,
                double x1max,
                int max2,
                double x2min,
                double x2max,
                std::function<double(double, double)> function2d,
                int romberg1,
                bool rational1,
                bool relative1,
                bool isLog1,
                int romberg2,
                bool rational2,
                bool relative2,
                bool isLog2,
                int rombergY,
                bool rationalY,
                bool relativeY,
                bool logSubst);

    //----------------------------------------------------------------------------//

    /*!
     * Constructor for the 1-dimensional functions if the array already exists.
     *
     * Constructs and initializes class for the 1-dimensional function with the
     * given array of x and f(x) values.
     *
     * \param   x    vector of x values.
     * \param   xmin vector of f(x) values
     * \param   romberg      Order of interpolation
     * \param   rational     interpolate with rational function
     * \param   relative     save error relative to the function value?
     * \return
     */
    Interpolant(std::vector<double> x, std::vector<double> y, int romberg, bool rational, bool relative);

    //----------------------------------------------------------------------------//

    /*!
     * Constructor for the 2-dimensional functions if the array already exists.
     *
     * Constructs and initializes class for the 2-dimensional function with the
     * given array of x1, x2 and f(x1, x2) values.
     *
     * \param   x1    vector of x values.
     * \param   x2    vector of x values.
     * \param   y vector of f(x) values
     * \param   romberg1      Order of interpolation
     * \param   rational1     interpolate with rational function
     * \param   relative1     save error relative to the function value?
     * \param   romberg2      Order of interpolation
     * \param   rational2     interpolate with rational function
     * \param   relative2     save error relative to the function value?
     * \return
     */
    Interpolant(std::vector<double> x1, std::vector<double> x2, std::vector< std::vector<double> > y, int romberg1, bool rational1, bool relative1 , int romberg2, bool rational2, bool relative2);

    //----------------------------------------------------------------------------//

    /*!
     * Constructor for the 2-dimensional functions if the array already exists, but x2 is asymmetric
     *
     * Constructs and initializes class for the 2-dimensional function with the
     * given array of x1, x2 and f(x1, x2) values.
     *
     * \param   x1    vector of x values.
     * \param   x2    vector of vectors with x2 values for every x1
     * \param   y vector of f(x) values
     * \param   romberg1      Order of interpolation
     * \param   rational1     interpolate with rational function
     * \param   relative1     save error relative to the function value?
     * \param   romberg2      Order of interpolation
     * \param   rational2     interpolate with rational function
     * \param   relative2     save error relative to the function value?
     * \return
     */
    Interpolant(std::vector<double> x1, std::vector< std::vector<double> > x2, std::vector< std::vector<double> > y, int romberg1, bool rational1, bool relative1 , int romberg2, bool rational2, bool relative2);

    //----------------------------------------------------------------------------//

    /**
     * Interpolates f(x) for 1d function
     *
     * \param    x
     * \return   interpolated value f(x)
     */

    double Interpolate(double x);

    //----------------------------------------------------------------------------//

    /**
     * Interpolates f(x) for 2d function
     *
     * \param    x1
     * \param    x2
     * \return   interpolated value f(x1,x2)
     */

    double Interpolate(double x1, double x2);

    //----------------------------------------------------------------------------//

    /**
     * Interpolates f(x) for 1d function if the arrays already exist.
     *
     * \param    x
     * \return   interpolated value f(x)
     */

    double InterpolateArray(double x);

    //----------------------------------------------------------------------------//

    /**
     * Interpolates f(x1,x2) for 2d function if the arrays already exist.
     *
     * \param    x1
     * \param    x2
     * \return   interpolated value f(x1,x2)
     */

    double InterpolateArray(double x1, double x2);

    //----------------------------------------------------------------------------//

    /**
     * Finds x: f(x)=y., 1d initialization required
     *
     * Finds x: f(x)=y. 1d initialization required.
     *
     * \param    y
     * \return   interpolated value x(y);
     */

    double FindLimit(double y);

    //----------------------------------------------------------------------------//

    /**
     * Finds x: f(a,x)=y.
     *
     * Finds x: f(a,x)=y, 2d initialization required
     *
     * \param    x1
     * \param    y
     * \return   interpolated value x(y);
     */

    double FindLimit(double x1, double y);

    //----------------------------------------------------------------------------//

    void swap(Interpolant& interpolant);

    //----------------------------------------------------------------------------//
    /**
     * Saves an interpolation table from file
     *
     * \param    Path/ofstream
     * \return   true if successfull
     */

    bool Save(std::string Path, bool binary_tables = false);
    bool Save(std::ofstream& out, bool binary_tables = false);

    //----------------------------------------------------------------------------//

    /**
     * Loads an interpolation table from file
     *
     * \param    Path/ifstream
     * \return   true if successfull
     */

    bool Load(std::string Path, bool binary_tables = false);
    bool Load(std::ifstream& in, bool binary_tables = false);

    //----------------------------------------------------------------------------//
    //----------------------------------------------------------------------------//
    //----------------------------------------------------------------------------//
    //----------------------------------------------------------------------------//
    // Getter
    /**
     * Getter for Interpolant object.
     *
     * \return   Interpolant object;
     */
    std::vector<Interpolant*> GetInterpolant() { return Interpolant_; }

    int GetRombergY() const { return rombergY_; }

    int GetRomberg() const { return romberg_; }

    std::vector<double> GetIX() const { return iX_; }

    std::vector<double> GetIY() const { return iY_; }

    std::vector<double> GetC() const { return c_; }

    std::vector<double> GetD() const { return d_; }

    int GetMax() const { return max_; }

    double GetXmin() const { return xmin_; }

    double GetXmax() const { return xmax_; }

    double GetStep() const { return step_; }

    bool GetRelative() const { return relative_; }

    bool GetRational() const { return rational_; }

    int GetRow() const { return row_; }

    int GetStarti() const { return starti_; }

    bool GetRationalY() const { return rationalY_; }

    bool GetRelativeY() const { return relativeY_; }

    bool GetSelf() const { return self_; }

    bool GetFlag() const { return flag_; }

    bool GetReverse() const { return reverse_; }

    bool GetIsLog() const { return isLog_; }

    bool GetLogSubst() const { return logSubst_; }

    double GetPrecision() const { return precision_; }

    double GetWorstX() const { return worstX_; }

    double GetWorstX2() const { return worstX2_; }

    double GetWorstY() const { return worstY_; }

    double GetPrecisionY() const { return precisionY_; }

    bool GetFast() const { return fast_; }

    double GetX_save() const { return x_save_; }

    double GetY_save() const { return y_save_; }

    //----------------------------------------------------------------------------//
    // Setter

    void SetRombergY(int rombergY);
    void SetRomberg(int romberg);
    void SetIX(const std::vector<double>& iX);
    void SetIY(const std::vector<double>& iY);
    void SetC(const std::vector<double>& c);
    void SetD(const std::vector<double>& d);
    void SetMax(int max);
    void SetXmin(double xmin);
    void SetXmax(double xmax);
    void SetStep(double step);
    void SetRelative(bool relative);
    void SetRational(bool rational);
    void SetRow(int row);
    void SetStarti(int starti);
    void SetRationalY(bool rationalY);
    void SetRelativeY(bool relativeY);
    void SetSelf(bool self);
    void SetFlag(bool flag);
    void SetReverse(bool reverse);
    void SetIsLog(bool isLog);
    void SetLogSubst(bool logSubst);
    void SetPrecision(double precision);
    void SetWorstX(double worstX);
    void SetWorstX2(double worstX2);
    void SetWorstY(double worstY);
    void SetPrecisionY(double precisionY);
    void SetFast(bool fast);
    void SetX_save(double x_save);
    void SetY_save(double y_save);
    /*!
     * Destructor
     */

    ~Interpolant();
};

} // namespace PROPOSAL
