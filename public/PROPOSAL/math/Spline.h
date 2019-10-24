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
