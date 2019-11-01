
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

#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include "PROPOSAL/math/Function.h"

namespace PROPOSAL {

struct spline_container {
   public:
    std::vector<double> coeff;
    std::pair<double, double> domain;

    friend std::ostream& operator<<(std::ostream& os, spline_container& s);
    friend std::istream& operator>>(std::istream& is, spline_container& s);
};

}  // namespace PROPOSAL

namespace PROPOSAL {

class Spline {
   public:
    Spline(std::vector<double>, std::vector<double>);
    Spline(std::vector<Polynom>, std::vector<double>);
    Spline(std::string spline_path, bool binary);
    Spline(const Spline&);
    virtual ~Spline() {};

    bool operator==(const Spline& spline) const;
    bool operator!=(const Spline& spline) const;

    virtual Spline* clone() const = 0;

    virtual bool save(std::string, bool);
    virtual double evaluate(double x);
    virtual void Derivative();
    virtual void Antiderivative(double c);

    std::pair<double, double> GetDomain();
    std::vector<Polynom> GetFunctions() { return splines_; }
    std::vector<double> GetSubintervalls() { return subintervall_; };
    unsigned int GetNumberIntervalls() { return n_subintervalls_; };

    virtual std::vector<spline_container> GetSplineContainer() const;

    friend std::ostream& operator<<(std::ostream& os, const Spline& s);
    friend std::istream& operator>>(std::istream& is, Spline& s);

   protected:
    virtual void calculate_splines(std::vector<double> x,
                                   std::vector<double> y) = 0;

    std::vector<Polynom> splines_;
    std::vector<double> subintervall_;
    unsigned int n_subintervalls_;
    std::vector<double> x_;
    std::vector<double> y_;
};

}  // namespace PROPOSAL

//---------------------------------------------------------------------------//
//------------------------------- Linear Spline -----------------------------//
//---------------------------------------------------------------------------//

namespace PROPOSAL {
class Linear_Spline : public Spline {
   public:
    Linear_Spline(std::vector<double>, std::vector<double>);
    Linear_Spline(std::vector<Polynom>, std::vector<double>);
    Linear_Spline(std::string spline_path, bool binary);
    Linear_Spline(const Spline&);

    Spline* clone() const override { return new Linear_Spline(*this); };

   protected:
    void calculate_splines(std::vector<double> x,
                           std::vector<double> y) override;
};

}  // namespace PROPOSAL

//----------------------------------------------------------------------------//
//-------------------------------- Cubic Spline
//------------------------------//
//----------------------------------------------------------------------------//

namespace PROPOSAL {

class Cubic_Spline : public Spline {
   public:
    Cubic_Spline(std::vector<double>, std::vector<double>);
    Cubic_Spline(std::vector<Polynom>, std::vector<double>);
    Cubic_Spline(std::string spline_path, bool binary);
    Cubic_Spline(const Spline&);

    Spline* clone() const override { return new Cubic_Spline(*this); };

   protected:
    void calculate_splines(std::vector<double> x,
                           std::vector<double> y) override;
};

}  // namespace PROPOSAL
